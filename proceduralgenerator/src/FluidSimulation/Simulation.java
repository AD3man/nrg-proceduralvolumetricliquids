package FluidSimulation;

/**
 * From: https://www.mikeash.com/pyblog/fluid-simulation-for-dummies.html
 */
public class Simulation {

    private FluidCube fluidCube;

    public Simulation(int size, double[][][] initialDensity, double viscosity, double diffusion, double dt) {
        fluidCube = new FluidCube(size,diffusion, viscosity, dt);
        for (int x = 0; x < size; x++) {
            for (int y = 0; y < size; y++) {
                for (int z = 0; z < size; z++) {
                    fluidCube.addDensity(x,y,z,initialDensity[x][y][z]);
                }
            }
        }

    }

    private void set_bnd(int b, double[][][] x, int N) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[i][j][0] = b == 3 ? -x[i][j][1] : x[i][j][1];
                x[i][j][N - 1] = b == 3 ? -x[i][j][N - 2] : x[i][j][N - 2];
            }
        }
        for (int k = 1; k < N - 1; k++) {
            for (int i = 1; i < N - 1; i++) {
                x[i][0][k] = b == 2 ? -x[i][1][k] : x[i][1][k];
                x[i][N - 1][k] = b == 2 ? -x[i][N - 2][k] : x[i][N - 2][k];
            }
        }
        for (int k = 1; k < N - 1; k++) {
            for (int j = 1; j < N - 1; j++) {
                x[0][j][k] = b == 1 ? -x[1][j][k] : x[1][j][k];
                x[N - 1][j][k] = b == 1 ? -x[N - 2][j][k] : x[N - 2][j][k];
            }
        }

        x[0][0][0] = 0.33f * (
                x[1][0][0]
                        + x[0][1][0]
                        + x[0][0][1]
        );
        x[0][N - 1][0] = 0.33f * (
                x[1][N - 1][0]
                        + x[0][N - 2][0]
                        + x[0][N - 1][1]
        );
        x[0][0][N - 1] = 0.33f * (
                x[1][0][N - 1]
                        + x[0][1][N - 1]
                        + x[0][0][N]
        );
        x[0][N - 1][N - 1] = 0.33f * (
                x[1][N - 1][N - 1]
                        + x[0][N - 2][N - 1]
                        + x[0][N - 1][N - 2]
        );
        x[N - 1][0][0] = 0.33f * (
                x[N - 2][0][0]
                        + x[N - 1][1][0]
                        + x[N - 1][0][1]
        );
        x[N - 1][N - 1][0] = 0.33f * (
                x[N - 2][N - 1][0]
                        + x[N - 1][N - 2][0]
                        + x[N - 1][N - 1][1]
        );
        x[N - 1][0][N - 1] = 0.33f * (
                x[N - 2][0][N - 1]
                        + x[N - 1][1][N - 1]
                        + x[N - 1][0][N - 2]
        );
        x[N - 1][N - 1][N - 1] = 0.33f * (
                x[N - 2][N - 1][N - 1]
                        + x[N - 1][N - 2][N - 1]
                        + x[N - 1][N - 1][N - 2]
        );
    }


    private void lin_solve(int b, double[][][] x, double[][][] x0, double a, double c, int iter, int N) {
        double cRecip = 1.0 / c;
        for (int k = 0; k < iter; k++) {
            for (int m = 1; m < N - 1; m++) {
                for (int j = 1; j < N - 1; j++) {
                    for (int i = 1; i < N - 1; i++) {
                        x[i][j][m] =
                                (x0[i][j][m]
                                        + a * (x[i + 1][j][m]
                                        + x[i - 1][j][m]
                                        + x[i][j + 1][m]
                                        + x[i][j - 1][m]
                                        + x[i][j][m + 1]
                                        + x[i][j][m - 1]
                                )) * cRecip;
                    }
                }
            }
            set_bnd(b, x, N);
        }
    }

    private void diffuse(int b, double[][][] x, double[][][] x0, double diffusion, double dt, int iter, int N) {
        double a = dt * diffusion * (N - 2) * (N - 2);
        lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
    }

    private void project(double[][][] velocX, double[][][] velocY, double[][][] velocZ, double[][][] p, double[][][] div, int iter, int N) {
        for (int k = 1; k < N - 1; k++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    div[i][j][k] = -0.5f * (
                            velocX[i + 1][j][k]
                                    - velocX[i - 1][j][k]
                                    + velocY[i][j + 1][k]
                                    - velocY[i][j - 1][k]
                                    + velocZ[i][j][k + 1]
                                    - velocZ[i][j][k - 1]
                    ) / N;
                    p[i][j][k] = 0;
                }
            }
        }
        set_bnd(0, div, N);
        set_bnd(0, p, N);
        lin_solve(0, p, div, 1, 6, iter, N);

        for (int k = 1; k < N - 1; k++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    velocX[i][j][k] -= 0.5f * (p[i + 1][j][k]
                            - p[i - 1][j][k]) * N;
                    velocY[i][j][k] -= 0.5f * (p[i][j + 1][k]
                            - p[i][j - 1][k]) * N;
                    velocZ[i][j][k] -= 0.5f * (p[i][j][k + 1]
                            - p[i][j][k - 1]) * N;
                }
            }
        }
        set_bnd(1, velocX, N);
        set_bnd(2, velocY, N);
        set_bnd(3, velocZ, N);
    }

    private void advect(int b, double[][][] d, double[][][] d0, double[][][] velocX, double[][][] velocY, double[][][] velocZ, double dt, int N) {
        double i0, i1, j0, j1, k0, k1;

        double dtx = dt * (N - 2);
        double dty = dt * (N - 2);
        double dtz = dt * (N - 2);

        double s0, s1, t0, t1, u0, u1;
        double tmp1, tmp2, tmp3, x, y, z;

        double Nfloat = N;
        double ifloat, jfloat, kfloat;
        int i, j, k;

        for (k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
            for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
                for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                    tmp1 = dtx * velocX[i][j][k];
                    tmp2 = dty * velocY[i][j][k];
                    tmp3 = dtz * velocZ[i][j][k];
                    x = ifloat - tmp1;
                    y = jfloat - tmp2;
                    z = kfloat - tmp3;

                    if (x < 0.5f) x = 0.5f;
                    if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
                    i0 = Math.floor(x);
                    i1 = i0 + 1.0f;
                    if (y < 0.5f) y = 0.5f;
                    if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
                    j0 = Math.floor(y);
                    j1 = j0 + 1.0f;
                    if (z < 0.5f) z = 0.5f;
                    if (z > Nfloat + 0.5f) z = Nfloat + 0.5f;
                    k0 = Math.floor(z);
                    k1 = k0 + 1.0f;

                    s1 = x - i0;
                    s0 = 1.0f - s1;
                    t1 = y - j0;
                    t0 = 1.0f - t1;
                    u1 = z - k0;
                    u0 = 1.0f - u1;

                    int i0i = (int) i0;
                    int i1i = (int) i1;
                    int j0i = (int) j0;
                    int j1i = (int) j1;
                    int k0i = (int) k0;
                    int k1i = (int) k1;

                    d[i][j][k] =

                            s0 * (t0 * (u0 * d0[i0i][j0i][k0i]
                                    + u1 * d0[i0i][j0i][k1i])
                                    + (t1 * (u0 * d0[i0i][j1i][k0i]
                                    + u1 * d0[i0i][j1i][k1i])))
                                    + s1 * (t0 * (u0 * d0[i1i][j0i][k0i]
                                    + u1 * d0[i1i][j0i][k1i])
                                    + (t1 * (u0 * d0[i1i][j1i][k0i]
                                    + u1 * d0[i1i][j1i][k1i])));
                }
            }
        }
        set_bnd(b, d, N);
    }

    public void step() {
        FluidCube cube = this.fluidCube;
        int N = cube.size - 1;
        double visc = cube.viscosity;
        double diff = cube.diffussion;
        double dt = cube.deltaTime;
        double[][][] Vx = cube.veloctityX;
        double[][][] Vy = cube.veloctityY;
        double[][][] Vz = cube.veloctityZ;
        double[][][] Vx0 = cube.veloctityXScratch;
        double[][][] Vy0 = cube.veloctityYScratch;
        double[][][] Vz0 = cube.veloctityZScratch;
        double[][][] s = cube.scatch;
        double[][][] density = cube.density;

        diffuse(1, Vx0, Vx, visc, dt, 4, N);
        diffuse(2, Vy0, Vy, visc, dt, 4, N);
        diffuse(3, Vz0, Vz, visc, dt, 4, N);

        project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);

        advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
        advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
        advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);

        project(Vx, Vy, Vz, Vx0, Vy0, 4, N);

        diffuse(0, s, density, diff, dt, 4, N);
        advect(0, density, s, Vx, Vy, Vz, dt, N);
    }

    public double[][][] getState() {
        return fluidCube.getDensityStateCopy();
    }
}
