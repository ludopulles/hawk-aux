from sage.all import load
from sys import argv
load("attack_Zn.sage")

experiments = int(argv[1])
processes = int(argv[2])

n = int(argv[3])

if n < 120:
    float_type = "double"
elif n < 210:
    float_type = "ld"
else:
    float_type = "qd"

sigma_min = float(argv[4])
sigma_max = float(argv[5])
sigma_step = float(argv[6])

g6k = bool(argv[7])

print("n, \t,sigma \t beta,\t prev_sd")

sigma = sigma_min

while sigma < sigma_max:
    beta, prev_sd, _, _ = many_experiment(experiments, processes, n,
                                          sigma=sigma, float_type=float_type,
                                          g6k=g6k)
    print("%d,%.3f,%.3f,%.3f" % (n, sigma, beta, prev_sd))
    sigma += sigma_step
