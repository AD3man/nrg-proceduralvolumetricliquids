// #package js/main

// #include WebGL.js

class Volume {

constructor(gl, reader, options) {
    Object.assign(this, {
        ready: false
    }, options);

    this._gl = gl;
    this._reader = reader;

    this.meta       = null;
    this.modalities = null;
    this.blocks     = null;
    this._texture   = null;
}

readMetadata(handlers) {
    if (!this._reader) {
        return;
    }
    this.ready = false;
    this._reader.readMetadata({
        onData: data => {
            this.meta = data.meta;
            this.modalities = data.modalities;
            this.blocks = data.blocks;
            handlers.onData && handlers.onData();
        }
    });
}

readModality(modalityName, handlers) {
    if (!this._reader || !this.modalities) {
        return;
    }
    this.ready = false;
    const modality = this.modalities.find(modality => modality.name === modalityName);
    if (!modality) {
        return;
    }
    const dimensions = modality.dimensions;
    const components = modality.components;
    const blocks = this.blocks;

    const gl = this._gl;
    if (this._texture) {
        gl.deleteTexture(this._texture);
    }
    this._texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_3D, this._texture);

    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_R, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);

        // TODO: here read modality format & number of components, ...
        let format, internalFormat;
        // if (components === 2) {
        //     internalFormat = gl.RG8;
        //     format = gl.RG;
        // } else {
        //     internalFormat = gl.R8;
        //     format = gl.RED;
        // }
        internalFormat = gl.RGBA8;
        format = gl.RGBA;


        gl.texStorage3D(gl.TEXTURE_3D, 1, internalFormat, dimensions.width, dimensions.height, dimensions.depth);
        let remainingBlocks = modality.placements.length;
        modality.placements.forEach(placement => {
            this._reader.readBlock(placement.index, {
                onData: data => {
                    let dataUint = new Uint8Array(data);
                    let newData = [];
                    let i2 = 0;

                    let gradients = [];
                    for (let k = 0; k < dimensions.depth; k++) {
                        for (let j = 0; j < dimensions.height; j++) {
                            for (let i = 0; i < dimensions.width; i++) {
                                let gradientsA = this.computeGradientSobelDistinct(
                                    dataUint,
                                    {
                                        width: dimensions.width,
                                        height: dimensions.height,
                                        depth: dimensions.depth
                                    },
                                    {
                                        x: i,
                                        y: j,
                                        z: k
                                    }
                                );
                                gradients.push(gradientsA);
                            }
                        }
                    }
                    for (let i = 0; i < data.byteLength; i++) {
                        //     // data[i] = -1;
                        //     newData[i2] = data[i];
                        //     newData[i2+1] = 0;
                        //     i2+=2;


                        newData.push(dataUint[i]);
                        newData.push(gradients[i].dx);
                        newData.push(gradients[i].dy);
                        newData.push(gradients[i].dz);

                   }
                   let newDataUint = Uint8Array.from(newData);

                     // debugger;
                    const position = placement.position;
                    const block = blocks[placement.index];
                    const blockdim = block.dimensions;
                    gl.bindTexture(gl.TEXTURE_3D, this._texture);
                    gl.texSubImage3D(gl.TEXTURE_3D, 0,
                        position.x, position.y, position.z,
                        blockdim.width, blockdim.height, blockdim.depth,
                        format, gl.UNSIGNED_BYTE, newDataUint);
                    remainingBlocks--;
                    if (remainingBlocks === 0) {
                        this.ready = true;
                        handlers.onLoad && handlers.onLoad();
                    }
                }
            });
        });
    }

    computeGradientSobelDistinct(volumeData, volumeDimensions, idx) {
        let kernelX = [
            1, 0, -1,
            2, 0, -2,
            1, 0, -1,

            2, 0, -2,
            4, 0, -4,
            2, 0, -2,

            1, 0, -1,
            2, 0, -2,
            1, 0, -1
        ];
        let kernelY = [
            1, 2, 1,
            0, 0, 0,
            -1, -2, -1,

            2, 4, 2,
            0, 0, 0,
            -2, -4, -2,

            1, 2, 1,
            0, 0, 0,
            -1, -2, -1
        ];
        let kernelZ = [
            1, 2, 1,
            2, 4, 2,
            1, 2, 1,

            0, 0, 0,
            0, 0, 0,
            0, 0, 0,

            -1, -2, -1,
            -2, -4, -2,
            -1, -2, -2
        ];
        let values = [
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,

            0, 0, 0,
            0, 0, 0,
            0, 0, 0,

            0, 0, 0,
            0, 0, 0,
            0, 0, 0
        ];
        for (let k = 0; k < 3; k++) {
            for (let j = 0; j < 3; j++) {
                for (let i = 0; i < 3; i++) {
                    let kernelIdx = Math.max(0,this.index({
                        x: i,
                        y: j,
                        z: k
                    }, {
                        width: 3,
                        height: 3,
                        depth: 3
                    }));
                    let volumeIndex = Math.max(0, this.index({
                        x: idx.x + i - 1,
                        y: idx.y + j - 1,
                        z: idx.z + k - 1
                    }, volumeDimensions));
                    values[kernelIdx] = volumeData[volumeIndex] / 255;
                }
            }
        }

        let dx = 0;
        for (let i = 0; i < kernelX.length; i++) {
            dx += kernelX[i] * values[i];
        }
        dx /= kernelX.length;

        let dy = 0;
        for (let i = 0; i < kernelY.length; i++) {
            dy += kernelY[i] * values[i];
        }
        dy /= kernelY.length;

        let dz = 0;
        for (let i = 0; i < kernelZ.length; i++) {
            dz += kernelZ[i] * values[i];
        }
        dz /= kernelZ.length;

        return {
            dx: dx,
            dy: dy,
            dz: dz
        };
    }

    index(idx, dims, offset) {
        offset = offset || {
            x: 0,
            y: 0,
            z: 0
        };
        return (idx.x + offset.x)
            + (idx.y + offset.y) * dims.width
            + (idx.z + offset.z) * dims.width * dims.height;
    }

    getTexture() {
        if (this.ready) {
            return this._texture;
        } else {
            return null;
        }
    }

    setFilter(filter) {
        if (!this._texture) {
            return;
        }

        var gl = this._gl;
        filter = filter === 'linear' ? gl.LINEAR : gl.NEAREST;
        gl.bindTexture(gl.TEXTURE_3D, this._texture);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, filter);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, filter);
    }

}
