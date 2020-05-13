import FluidSimulation.Simulation;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.io.*;
import java.util.ArrayList;
import java.util.Random;

public class VPTWriterTest {

    public final int X = 128;
    public final int Y = 128;
    public final int Z = 128;

    public double[][][] prepareInitialConditions(double waterDensity, int seed) {
        double[][][] data = new double[X][Y][Z];
        Random r = seed > 0 ? new Random(seed) : new Random();
        for (int x = 0; x < X; x++) {
            for (int y = 0; y < Y; y++) {
                for (int z = 0; z < Z; z++) {
                    double a = ImprovedNoise.noise(x + r.nextDouble(), y + r.nextDouble(), z + r.nextDouble());
                    data[x][y][z] = waterDensity + (100 * a);
                }
            }
        }
        return data;
    }

    public double[][][] doSimulation(double[][][] initial, double viscosity, double diffusion, double dt) {
        Simulation s = new Simulation(X, initial, viscosity, diffusion, dt);
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
                    double s = raster.getSample(i,j,0);
                    double s2 = (s)/255 * ((double)(Z-1));

                   imgArr[i][j] = (int)Math.round(s2);
                    if( imgArr[i][j] > 0) {
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
            }
        }
        return reducedSeaLevel;
    }

    private double[][][] addGround(double[][][] simulated, int[][] heightmap, int meanZ, int dZ, double density, int seed) {
        double[][][] withGround = new double[X][Y][Z];
        Random r = seed > 0 ? new Random(seed) : new Random();
        for (int x = 0; x < X; x++) {
            for (int y = 0; y < Y; y++) {
                double a;
                int zMax;
                if (heightmap == null) {
                    a = ImprovedNoise.noise(x + r.nextDouble(), y + r.nextDouble(), meanZ);
                    zMax = (int) Math.round(meanZ + (a * dZ));
                } else {
                    a = ImprovedNoise.noise(x + r.nextDouble(), y + r.nextDouble(), meanZ + r.nextDouble());
                    zMax = Math.max(0,Math.min((int) Math.round(heightmap[x][y] + a), Z-1));
                }

                for (int z = 0; z < zMax; z++) {
                    withGround[x][y][z] = density;
                }
                for (int z = zMax; z < Z; z++) {
                    withGround[x][y][z] = simulated[x][y][z];
                }
            }
        }
        return withGround;
    }

    public void writeToFile(String path, double[][][] data, double maxDensity, double minDensity) {

        maxDensity = 1000 * maxDensity;
        minDensity = 1000 * minDensity;

        try (FileOutputStream fo = new FileOutputStream(path);) {
            ArrayList<Integer> noise = new ArrayList<>();
            Random r = new Random();
            for (int x = 0; x < X; x++) {
                for (int y = 0; y < Y; y++) {
                    for (int z = 0; z < Z; z++) {
                        double a = (((data[x][y][z]) - minDensity) / maxDensity) * 255;
                        if (a < 2) {
                            System.out.println("very small");
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
                    System.out.println("very small");
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
        double dt = 0.001;
        double viscosity = 1.5673;

        int[][] heightmap = vptw.getHeightmap("heightmaps/hmap3.jpg");
        String path = "resout/test5.bin";


        double[][][] initial = vptw.prepareInitialConditions(1000, 234);
        double[][][] simulated1 = vptw.doSimulation(initial, viscosity, diffusion, dt);
        double[][][] simulated2 = vptw.doSimulation(simulated1, viscosity, diffusion, dt);
        double[][][] simulated3 = vptw.doSimulation(simulated2, viscosity, diffusion, dt);
        double[][][] reducedSeaLevel = vptw.reduceSeaLevel(simulated2, 50);
//        double[][][] withGround = vptw.addGround(simulated3, null,10, 5, 2700, 235);
        double[][][] withGround = vptw.addGround(reducedSeaLevel, heightmap , 10, 0, 2700, 235);


        vptw.writeToFile(path, withGround, 2.7, 0.0012);
    }


}
