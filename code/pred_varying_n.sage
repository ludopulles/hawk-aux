from sage.all import load, log, sqrt
load("BKZ_simulator.sage")

print("n,\t beta,\t prev_sd")

ns = [i for i in range(50, 201)]
ns += [i for i in range(210, 451, 10)]
ns += [i for i in range(500, 851, 50)]

# for blocksizes and prev_sds of parameter sets
# ns = [512, 1024, 2048]

for n in ns:

    if n <= 250:
        beta, prev_sd = predict_beta_and_prev_sd(n, n*log(n)/2,
                                                 lift_union_bound=False,
                                                 number_targets=n, tours=1)
    else:
        # precision issues with non union bound methods, which
        # coincide with union bound estimates by this large blocksize
        beta, prev_sd = predict_beta_and_prev_sd(n, n*log(n)/2,
                                                 lift_union_bound=True,
                                                 number_targets=n, tours=1)

    # lattice was scaled up, so prev_sd need to be scaled back
    prev_sd /= sqrt(n)

    print("%d, \t%.3f,\t%.3f" % (n, beta, prev_sd))
