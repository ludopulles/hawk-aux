from proba_utils import centered_discrete_Gaussian_law, variance
from math import sqrt

# used to convert the Gaussian parameter sigma into a standard deviation
# tilde{sigma} for Fig 4 when sigma is below smoothing

for d in range(100, 190, 10):
    inp = open("../data/%d_varying_sigma.txt" % d, "rt")
    out = open("../data/%d_varying_stddev.txt" % d, "wt")
    out.write("n, \tstddev,  \tbeta,	\tprev_sd \n")
    for line in inp:
        try:
            exec("X =" + line)
            _, sigma, beta, prev_sd = X # noqa
        except NameError:
            # skip first line "n, stddev, beta, prev_sd"
            continue

        D = centered_discrete_Gaussian_law(sigma)
        stddev = sqrt(variance(D))
        out.write("%d, \t%.3f, \t%.2f, \t%.3f\n" % (d, stddev, beta, prev_sd))
