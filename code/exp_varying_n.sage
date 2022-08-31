from sage.all import load
from sys import argv

load("attack_Zn.sage")

experiments = int(argv[1])
processes = int(argv[2])

n_min = int(argv[3])
n_max = int(argv[4])
n_step = int(argv[5])

# 0 for unstructred, 1 for structured
structured = bool(int(argv[6]))

# 0 for FPyLLL, 1 for g6k
g6k = bool(int(argv[7]))

# 1 checks that both f, g have odd norm, as in Hawk
hawk = bool(int(argv[8]))

print("n,\t beta,\t prev_sd")

if structured:
    power = n_min
    while power < n_max:
        n = 2**power
        if 2*n <= 128:
            float_type = "double"
        else:
            float_type = "qd"

        beta, prev_sd, geo_prev_sd, _ = many_experiment_structured(experiments,
                                                                   processes,
                                                                   n,
                                                                   float_type=float_type,  # noqa
                                                                   g6k=g6k,
                                                                   hawk=hawk)

        data = "%d,%.3f,%.3f,%.3f\n" % (2*n, beta, prev_sd, geo_prev_sd)
        print("%d,\t %.3f,\t %.3f,\t %.3f" % (2*n, beta, prev_sd, geo_prev_sd))

        if g6k:
            filename = "structured-prev_sd-g-{nmin}-{nmax}-{nstep}-{experiments}".format(  # noqa
                        nmin=n_min, nmax=n_max, nstep=n_step,
                        experiments=experiments)
        else:
            filename = "structured-prev_sd-{nmin}-{nmax}-{nstep}-{experiments}".format(  # noqa
                        nmin=n_min, nmax=n_max, nstep=n_step,
                        experiments=experiments)

        with open("../data/" + filename, "a") as out:
            out.write(data)
            out.close()

        power += n_step
else:
    n = n_min
    while n < n_max:
        if n < 120:
            float_type = "double"
        elif n < 210:
            float_type = "ld"
        else:
            float_type = "qd"
        beta, prev_sd, geo_prev_sd, _ = many_experiment(experiments, processes,
                                                        n,
                                                        float_type=float_type,
                                                        g6k=g6k)

        data = "%d,%.3f,%.3f,%.3f\n" % (n, beta, prev_sd, geo_prev_sd)
        print("%d,\t %.3f,\t %.3f,\t %.3f" % (n, beta, prev_sd, geo_prev_sd))

        if g6k:
            filename = "prev_sd-g-{nmin}-{nmax}-{nstep}-{experiments}".format(
                        nmin=n_min, nmax=n_max, nstep=n_step,
                        experiments=experiments)
        else:
            filename = "prev_sd-{nmin}-{nmax}-{nstep}-{experiments}".format(
                        nmin=n_min, nmax=n_max, nstep=n_step,
                        experiments=experiments)

        with open("../data/" + filename, "a") as out:
            out.write(data)
            out.close()

        n += n_step
