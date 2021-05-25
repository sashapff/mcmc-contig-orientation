from scipy.integrate import quad
from scipy.stats import expon, uniform, norm, gamma
from tqdm import tqdm
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def uniform_naive_distribution(n_reads=10000, scale=100, D=50):
    distances = np.zeros(n_reads, dtype=np.int64)

    for i in tqdm(range(n_reads)):
        distance = int(uniform.rvs(scale=scale))
        distances[i] = distance

    distances_D = distances[distances > D]

    P = lambda x: np.where(x < 0, 0, np.where(x < scale, 1 / (2 * distances.mean()), 0))

    quad_value, _ = quad(lambda x: P(x), 0, np.inf)
    print(f'quad_value: {quad_value}')

    sns.distplot(distances, norm_hist=True, kde=False)
    plt.title(f'Real distances distribution')
    x_range = range(distances.max())
    plt.plot(x_range, [P(x) for x in x_range])
    plt.show()

    sns.distplot(distances_D, norm_hist=True, kde=False)
    plt.title(f'Distances distribution, d > {D}')
    plt.show()

    log_likelihoods = []

    d_range = range(0, 2 * D)
    for d in d_range:
        log_likelihood = 0
        quad_value, _ = quad(lambda x: P(x), d, np.inf)
        P_d = lambda x: np.where(x <= d, 0, P(x) / quad_value)
        for distance in distances_D - D + d:
            if distance > d:
                log_likelihood += np.log(P_d(distance))
            else:
                log_likelihood = -np.inf
        log_likelihoods.append(log_likelihood)

    plt.plot(d_range, log_likelihoods)
    plt.title('Log likelihood')
    plt.show()

def expon_naive_distribution(n_reads=10000, scale=100, D=50):
    distances = np.zeros(n_reads, dtype=np.int64)

    for i in tqdm(range(n_reads)):
        distance = int(expon.rvs(scale=scale))
        distances[i] = distance

    distances_D = distances[distances > D]

    Lambda = 1 / distances.mean()
    P = lambda x: Lambda * np.exp(-Lambda * x)

    quad_value, _ = quad(lambda x: P(x), 0, np.inf)
    print(f'quad_value: {quad_value}')

    log_likelihoods = []

    d_range = range(0, 2 * D)
    for d in d_range:
        log_likelihood = 0
        quad_value, _ = quad(lambda x: P(x), d, np.inf)
        P_d = lambda x: np.where(x <= d, 0, P(x) / quad_value)
        for distance in distances_D:
            if distance > d:
                log_likelihood += np.log(P_d(distance))
            else:
                log_likelihood = -np.inf
        log_likelihoods.append(log_likelihood)

    plt.plot(d_range, log_likelihoods)
    plt.title('Log likelihood for naive estimation')
    plt.xlabel('d, distance estimate')
    plt.ylabel('log likelihood for d')
    plt.show()

def norm_naive_distribution(n_reads=10000, scale=100, D=50):
    distances = np.zeros(n_reads, dtype=np.int64)

    for i in tqdm(range(n_reads)):
        distance = int(norm.rvs(scale=scale))
        distances[i] = distance

    distances_D = distances[distances > D]

    P = lambda x: norm.pdf(x, distances.mean(), distances.std())

    quad_value, _ = quad(lambda x: P(x), 0, np.inf)
    print(f'quad_value: {quad_value}')

    sns.distplot(distances, norm_hist=True, kde=False)
    plt.title(f'Real distances distribution')
    x_range = range(distances.min(), distances.max())
    plt.plot(x_range, [P(x) for x in x_range])
    plt.show()

    sns.distplot(distances_D, norm_hist=True, kde=False)
    plt.title(f'Distances distribution, d > {D}')
    plt.show()

    log_likelihoods = []

    d_range = range(-2 * D, D)
    for d in d_range:
        log_likelihood = 0
        quad_value, _ = quad(lambda x: P(x), d, np.inf)
        P_d = lambda x: np.where(x <= d, 0, P(x) / quad_value)
        for distance in distances_D - D + d:
            if distance > d:
                log_likelihood += np.log(P_d(distance))
            else:
                log_likelihood = -np.inf
        log_likelihoods.append(log_likelihood)

    plt.plot(d_range, log_likelihoods)
    plt.title('Log likelihood')
    plt.show()

def gamma_naive_distribution(n_reads=10000, scale=100, D=50):
        distances = np.zeros(n_reads, dtype=np.int64)

        for i in tqdm(range(n_reads)):
            distance = int(gamma.rvs(a=2, scale=scale))
            distances[i] = distance

        distances_D = []
        for d in distances:
            if d > D:
                distances_D.append(d)
            else:
                distances_D.append(1500)
        distances_D = np.array(distances_D)

        P = lambda x: gamma.pdf(x, a=2, scale=scale)

        quad_value, _ = quad(lambda x: P(x), -np.inf, np.inf)
        print(f'quad_value: {quad_value}')

        sns.distplot(distances_D, norm_hist=True, kde=False)
        plt.title(f'Observed distances distribution')
        x_range = range(distances.min(), distances.max())
        plt.plot(x_range, [P(x) for x in x_range], label='real probability density')
        plt.xlabel('d, distance')
        plt.axvline(x=50, label='correct distance', color='red')
        plt.legend()
        plt.xlim((distances.min(), distances.max()))
        plt.ylabel('p(d), probability density')
        plt.show()

        sns.distplot(distances, norm_hist=True, kde=False)
        plt.title(f'Real distances distribution')
        x_range = range(distances.max())
        plt.plot(x_range, [P(x) for x in x_range])
        plt.xlabel('d, distance')
        plt.ylabel('p(d), probability density')
        plt.xlim((distances.min(), distances.max()))
        plt.show()

        quad_value, _ = quad(lambda x: P(x), D, np.inf)
        P_D = lambda x: np.where(x <= D, 0, P(x) / quad_value)

        log_likelihoods = []

        d_range = range(0, 2 * D)
        for d in d_range:
            log_likelihood = 0
            quad_value, _ = quad(lambda x: P(x), d, np.inf)
            P_d = lambda x: np.where(x <= d, 0, P(x) / quad_value)
            for distance in distances_D - D + d:
                if distance > d:
                    log_likelihood += np.log(P_d(distance))
                else:
                    log_likelihood = -np.inf
            log_likelihoods.append(log_likelihood)
            if d % 25 == 0:
                plt.plot(x_range, [P_d(x) for x in x_range], label=f'd={d}')

        plt.legend()
        plt.show()

        plt.plot(d_range, log_likelihoods)
        plt.axvline(x=50, label='correct distance', color='red')
        plt.title('Log likelihood')
        plt.legend()
        plt.ylabel('log_likelihood(d)')
        plt.xlabel('d, distance estimate')
        plt.show()

def expon_real_distribution(contig_length=10_000, n_reads=10000, scale=1000, D=500):
    gap_start = contig_length // 2 - D // 2
    lefts = np.zeros(n_reads, dtype=np.int64)
    rights = np.zeros(n_reads, dtype=np.int64)
    distances = np.zeros(n_reads, dtype=np.int64)

    for i in tqdm(range(n_reads)):
        left = np.random.randint(0, contig_length)
        right = left + expon.rvs(scale=scale)
        while right >= contig_length:
            right = left + expon.rvs(scale=scale)
        lefts[i] = left
        rights[i] = right
        distances[i] = right - left

    Lambda = 1 / distances.mean()
    P = lambda x: np.where(x >= 0, Lambda * np.exp(-Lambda * x), 0)

    quad_value, _ = quad(lambda x: P(x), -np.inf, np.inf)
    print(f'quad_value: {quad_value}')

    distances_D = []
    for i in range(n_reads):
        if lefts[i] < gap_start and rights[i] > gap_start + D:
            distances_D.append(distances[i])
        else:
            distances_D.append(np.random.randint(contig_length, contig_length + 100))

    distances_D = np.array(distances_D)

    sns.distplot(distances, norm_hist=True, kde=False)
    plt.title(f'Real distances distribution')
    x_range = range(distances.min(), distances.max())
    plt.plot(x_range, [P(x) for x in x_range])
    plt.xlabel('d, distance')
    plt.ylabel('p(d), probability density')
    plt.show()

    sns.distplot(distances_D, norm_hist=True, kde=False)
    plt.title(f'Real distances distribution')
    x_range = range(distances_D.min(), distances_D.max())
    plt.plot(x_range, [P(x) for x in x_range], label='real probability density')
    plt.xlabel('d, distance')
    plt.ylabel('p(d), probability density')
    plt.show()

    sns.distplot(distances_D, norm_hist=True, kde=False)
    plt.title(f'Distances distribution, d > {D}')
    plt.show()

    log_likelihoods = []

    d_range = range(0, 2 * D)
    for d in d_range:
        log_likelihood = 0
        quad_value, _ = quad(lambda x: P(x), d, np.inf)
        P_d = lambda x: np.where(x <= d, 0, P(x) / quad_value)
        for distance in (distances_D - D + d):
            if distance > d:
                log_likelihood += np.log(P_d(distance))
            else:
                log_likelihood = -np.inf
        log_likelihoods.append(log_likelihood)

    plt.plot(d_range, log_likelihoods)
    plt.title('Log likelihood')
    plt.show()

def gamma_real_distribution(contig_length=10_000, n_reads=10000, scale=1000, D=500):
    gap_start = contig_length // 2 - D // 2
    lefts = []
    rights = []
    distances = []

    for i in tqdm(range(n_reads)):
        left = np.random.randint(0, contig_length)
        right = left + int(gamma.rvs(a=2, scale=scale))
        if right < contig_length:
            lefts.append(left)
            rights.append(right)
            distances.append(right - left)

    lefts = np.array(lefts)
    rights = np.array(rights)
    distances = np.array(distances)

    P = lambda x: gamma.pdf(x, a=2, scale=scale)

    quad_value, _ = quad(lambda x: P(x), -np.inf, np.inf)
    print(f'quad_value: {quad_value}')

    distances_D = []
    for i in range(len(lefts)):
        if lefts[i] < gap_start and rights[i] > gap_start + D:
            distances_D.append(distances[i])

    distances_D = np.array(distances_D)

    sns.distplot(distances, norm_hist=True, kde=False)
    plt.title(f'Real distances distribution')
    x_range = range(distances.min(), distances.max())
    plt.plot(x_range, [P(x) for x in x_range])
    plt.xlabel('d, distance')
    plt.ylabel('p(d), probability density')
    plt.show()

    sns.distplot(distances_D, norm_hist=True, kde=False)
    plt.title(f'Real distances distribution')
    x_range = range(distances_D.min(), distances_D.max())
    plt.plot(x_range, [P(x) for x in x_range], label='real probability density')
    plt.xlabel('d, distance')
    plt.ylabel('p(d), probability density')
    plt.show()


    log_likelihoods = []

    d_range = range(0, 2 * D, 25)
    for d in tqdm(d_range):
        log_likelihood = 0
        quad_value, _ = quad(lambda x: P(x), d, np.inf)
        P_d = lambda x: np.where(x <= d, 0, P(x) / quad_value)
        for distance in (distances_D - D + d):
            if distance > d:
                log_likelihood += np.log(P_d(distance))
            else:
                log_likelihood = -np.inf
        log_likelihoods.append(log_likelihood)

        sns.distplot(distances_D - D + d, norm_hist=True, kde=False)
        plt.title(f'Distances distribution, d={d}')

        plt.plot(x_range, [P_d(x) for x in x_range], label=f'd={d}')
        plt.xlabel('d, distance')
        plt.ylabel('p(d), probability density')
        plt.legend()
        plt.show()

    plt.plot(d_range, log_likelihoods)
    plt.title('Log likelihood')
    plt.show()

if __name__ == "__main__":
    gamma_naive_distribution()


