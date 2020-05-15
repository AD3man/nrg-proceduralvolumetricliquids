import FluidSimulation.Simulation;
import com.flowpowered.noise.Noise;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.io.*;
import java.util.ArrayList;
import java.util.Random;
import java.util.function.DoubleToIntFunction;

public class VPTWriterTest {

    public final int X = 128;
    public final int Y = 128;
    public final int Z = 128;

    public double[][][] prepareInitialConditions(double waterDensity, double zToRealCoef, double maxDensityDiff, double offsetDelta, int seed) {
        double[][][] data = new double[X][Y][Z];
        double g = 9.81;
        Random r = new Random();
        int offsetRanX = r.nextInt();
        int offsetRanY = r.nextInt();
        int offsetRanZ = r.nextInt();
        double xOffset = offsetRanX + 0.1;
        for (int x = 0; x < X; x++) {
            double yOffset = offsetRanY +0.1;
            for (int y = 0; y < Y; y++) {
                double zOffset = offsetRanZ +  0.1;
                for (int z = 0; z < Z; z++) {
                    double a = ImprovedNoise.noise(xOffset, yOffset, zOffset);
                    double h = Z - z;
                    data[x][y][z] = (waterDensity * g * (h * zToRealCoef)) + ((maxDensityDiff) * a);
                    zOffset += offsetDelta;
                }
                yOffset += offsetDelta;
            }
            xOffset += offsetDelta;
        }
        return data;
    }

    public double[][][] doSimulation(double[][][] initial, double viscosity, double diffusion, double dt, int seed) {
        Simulation s = new Simulation(X, initial, viscosity, diffusion, dt, seed);
        s.step();
        return s.getState();
    }

    private int[][] getHeightmap(String path) {
        File file = new File(path);
        try {
            BufferedImage img = ImageIO.read(file);
            int[][] imgArr = new int[X][Y];
            Raster raster = img.getData();
            for (int i = 0; i < X; i++) {
                for (int j = 0; j < Y; j++) {
                    double s = raster.getSample(i, j, 0);
                    double s2 = (s) / 255 * ((double) (Z - 1));

                    imgArr[i][j] = (int) Math.round(s2);
                    if (imgArr[i][j] > 0) {
                        int a = 1;
                    }
                }
            }
            return imgArr;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;

    }

    private double[][][] reduceSeaLevel(double[][][] simulated, int reduceSize) {
        double[][][] reducedSeaLevel = new double[X][Y][Z];
        int newZ = Z - reduceSize;
        for (int x = 0; x < X; x++) {
            for (int y = 0; y < Y; y++) {
                for (int z = 0; z < newZ; z++) {
                    reducedSeaLevel[x][y][z] = simulated[x][y][z];
                }
                for (int z = newZ; z < Z; z++) {
                    reducedSeaLevel[x][y][z] = 1.2;
                }

            }
        }
        return reducedSeaLevel;
    }

    private double[][][] reduceSeaLevelAddWaves(double[][][] simulated, int reduceSize, int waveAmpl, double offsetDelta) {
        double[][][] reducedSeaLevel = new double[X][Y][Z];
        int newZ = (Z - reduceSize);
        double newZd = newZ + 0.01;
//       TODO potisni pressure dol.
        Random r = new Random();
        int offsetRanX = r.nextInt();
        int offsetRanY = r.nextInt();
        double xOffset = offsetRanX + 0.1;
        for (int x = 0; x < X; x++) {
            double yOffset = offsetRanY + 0.1;
            for (int y = 0; y < Y; y++) {
                double a = ImprovedNoise.noise(xOffset, yOffset, newZd);
                int zMax = (int) Math.round(newZ + (waveAmpl * a));
                for (int z = 0; z < zMax; z++) {
                    if(z == zMax-1) {
                        double qwq = 1.0;
                    }
                    int sourceZ = z + (reduceSize + (newZ - zMax));
                    double value =  simulated[x][y][sourceZ];
                    reducedSeaLevel[x][y][z] = value;
                }
                for (int z = zMax; z < Z; z++) {
                    reducedSeaLevel[x][y][z] = 1.2;
                }
                yOffset += offsetDelta;
            }
            xOffset += offsetDelta;
        }
        return reducedSeaLevel;
    }

    private double[][][] addGround(double[][][] simulated, int[][] heightmap, int meanZ, int dZ, double offsetDelta, double density) {
        double[][][] withGround = new double[X][Y][Z];
        double meanZd = meanZ + 0.01;
        Random r = new Random();
        int offsetRanX = r.nextInt();
        int offsetRanY = r.nextInt();
        double xOffset =offsetRanX + 0.1;
        double g = 9.81;
        for (int x = 0; x < X; x++) {
            double yOffset = offsetRanY + 0.1;
            for (int y = 0; y < Y; y++) {
                double a;
                int zMax;
                if (heightmap == null) {
                    a = ImprovedNoise.noise(xOffset, yOffset, meanZd);
                    zMax = (int) Math.round(meanZ + (a * dZ));
                } else {
                    a = ImprovedNoise.noise(xOffset, yOffset, meanZd);
                    zMax = Math.max(0, Math.min((int) Math.round(heightmap[x][y] + (a * dZ)), Z - 1));
                }

                for (int z = 0; z < zMax; z++) {

                    withGround[x][y][z] = density;
                }
                for (int z = zMax; z < Z; z++) {
                    withGround[x][y][z] = simulated[x][y][z];
                }
                yOffset += offsetDelta;

            }
            xOffset += offsetDelta;
        }
        return withGround;
    }

    public void writeToFile(String path, String extension, double[][][] data, double maxDensity, double minDensity) {

        maxDensity = 1000 * maxDensity;
        minDensity = 1000 * minDensity;

        try (FileOutputStream fo = new FileOutputStream(path + extension);) {
            ArrayList<Integer> noise = new ArrayList<>();
            Random r = new Random();
            for (int x = 0; x < X; x++) {
                for (int y = 0; y < Y; y++) {
                    for (int z = 0; z < Z; z++) {
                        double a = (((data[x][y][z]) - minDensity) / maxDensity) * 255;
                        if (a < 2) {
//                            System.out.println("very small");
                        }
                        int bytev = (int) Math.round(a);

//                        noise[x][y][z] = bytev;
                        noise.add(bytev);
                    }
                }
            }
            int[] out = new int[X * Y * Z];
            for (int i = 0; i < X * Y * Z; i++) {
                out[i] = noise.get(i);
                if (out[i] < 2) {
//                    System.out.println("very small");
                }
                fo.write(out[i]);
            }

//            fo.write(out);
            System.out.println("done");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public static void main(String[] args) {
        VPTWriterTest vptw = new VPTWriterTest();
        double diffusion = 0.001;
//        double dt = 0.001;
//        double dt = 0.01;
        double dt = 0.01;
        double viscosity = 1.5673;

        int[][] heightmap = null; //vptw.getHeightmap("heightmaps/hmap3.jpg");
        String path = "resout/test012_dt0_01_velocity5_init";


        double[][][] initial = vptw.prepareInitialConditions(1000, 1.0 / 1000.0, 0.5, 0.01, 234);

        int startingseed = 235;
        double[][][] simulated = vptw.doSimulation(initial, viscosity, diffusion, dt, startingseed++);
        for (int step = 1; step < 10; step++) {
            simulated = vptw.doSimulation(simulated, viscosity, diffusion, dt, startingseed++);
            if (step % 2 == 0) {
                vptw.writeToFile(path + "_simstep_" + step, ".bin", simulated, 2.7, 0.0012);
            }
        }


        double maxP = Double.MIN_VALUE;
        double minP = Double.MAX_VALUE;
        for (int x = 0; x < vptw.X; x++) {
            for (int y = 0; y < vptw.Y; y++) {
                for (int z = 0; z < vptw.Z; z++) {
                    double value = simulated[x][y][z];
                    if (value > maxP) {
                        maxP = value;
                    }
                    if (value < minP) {
                        minP = value;
                    }
                }
            }
        }
        System.out.println(String.format("Simulated pressure min: %f; max: %f ", minP, maxP));


//        double[][][] reducedSeaLevel = vptw.reduceSeaLevel(simulated, 50);
//        vptw.writeToFile(path + "_reduce.bin", reducedSeaLevel, 2.7, 0.0012);
        path += "_reduce_waves_lowerdelta";
        double[][][] reducedSeaLevelAndWave = vptw.reduceSeaLevelAddWaves(simulated, 50, 10, 0.05);
        vptw.writeToFile(path, ".bin", reducedSeaLevelAndWave, 2.7, 0.0012);

//        double[][][] withGround = vptw.addGround(reducedSeaLevel, heightmap, 10, 0, 2700, 235);
//        double[][][] withGround = vptw.addGround(reducedSeaLevelAndWave, heightmap, 20, 5, 0.05, 2700);
        double[][][] withGround = vptw.addGround(reducedSeaLevelAndWave, heightmap, 20, 5, 0.05, maxP*2);


        System.out.println(String.format("Perlin noise min: %f; max: %f ", ImprovedNoise.min, ImprovedNoise.max));
//        vptw.writeToFile(path + "_ground", ".bin", withGround, 2.7, 0.0012);
        vptw.writeToFile(path + "_ground", ".bin", withGround, (maxP*2)/1000, 1.0/1000);
    }


}
