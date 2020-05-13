package FluidSimulation;
/**
 * From: https://www.mikeash.com/pyblog/fluid-simulation-for-dummies.html
 */
public class FluidCube {
    int size;
    double deltaTime;
    double diffussion;
    double viscosity;

    double[][][] scatch;
    double[][][] density;

    double[][][] veloctityX;
    double[][][] veloctityY;
    double[][][] veloctityZ;

    double[][][] veloctityXScratch;
    double[][][] veloctityYScratch;
    double[][][] veloctityZScratch;

    public FluidCube(int size, double diffussion, double viscosity, double dt) {
        this.size = size;
        this.deltaTime = dt;
        this.diffussion = diffussion;
        this.viscosity = viscosity;

        this.scatch = new double[size][size][size];
        this.density = new double[size][size][size];
        this.veloctityX = new double[size][size][size];
        this.veloctityY = new double[size][size][size];
        this.veloctityZ = new double[size][size][size];
        this.veloctityXScratch = new double[size][size][size];
        this.veloctityYScratch = new double[size][size][size];
        this.veloctityZScratch = new double[size][size][size];
    }

    public void addDensity(int x, int y, int z, double amount) {
        this.density[x][y][z] += amount;
    }

    public void addVelocity(int x, int y, int z, double amountX, double amountY, double amountZ) {
        this.veloctityX[x][y][z] += amountX;
        this.veloctityY[x][y][z] += amountY;
        this.veloctityZ[x][y][z] += amountZ;
    }

    public double[][][] getDensityStateCopy() {
        double[][][] state = new double[size][size][size];
        for(int x = 0; x < size; x++) {
            for(int y = 0; y < size; y++) {
                for(int z = 0; z < size; z++) {
                    state[x][y][z] = density[x][y][z];
                }
            }
        }
        return state;
    }
}
