from scipy.stats import expon

from orientation.tools import simulation

if __name__ == "__main__":
    simulation(max_len=1000, n_reads=100, n_contigs=5, p_distr=expon, output_path='../data_sim')