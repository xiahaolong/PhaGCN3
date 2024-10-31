import os
import numpy as np
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import gc

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--outpath', type=str, default="result")
args = parser.parse_args()

def hypergeom_sf(x, N, K, n):
    p = 0
    for k in range(int(x) + 1):
        if k < 0 or k > K or (n - k) < 0 or (N - K) < (n - k):
            continue  # 跳过无效的 k 值
        p += (math.comb(int(K), int(k)) * math.comb(int(N - K), int(n - k))) / math.comb(int(N), int(n))
    return 1 - p

def calculate_significance(A, B, commons_pc, pcs, number_of_pc, logT, max_sig):
    a, b = sorted([number_of_pc[A], number_of_pc[B]])
    commons_value = commons_pc[A, B] - 1
    
    # 确保 commons_value 非负
    if commons_value < 0:
        return A, B, 0  # 返回0表示不显著
    
    pval = hypergeom_sf(commons_value, pcs, a, b)
    sig = min(max_sig, np.nan_to_num(-np.log10(pval) - logT))
    return A, B, sig

def create_network(matrix_path, singletons_path, thres=1, max_sig=1000, chunk_size=100):
    # 使用内存映射加载数据
    matrix = np.load(matrix_path, mmap_mode='r')
    singletons = np.load(singletons_path, mmap_mode='r')

    contigs = matrix.shape[0]
    pcs = matrix.shape[1] + singletons.sum()
    T = 0.5 * contigs * (contigs - 1)
    logT = np.log10(T)
    number_of_pc = matrix.sum(axis=1) + singletons
    commons_pc = np.dot(matrix, matrix.T)

    S = np.zeros((contigs, contigs), dtype=np.float32)

    # 分块处理
    for start in range(0, contigs, chunk_size):
        end = min(start + chunk_size, contigs)
        with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
            futures = []
            for A in range(start, end):
                for B in range(A + 1, contigs):
                    if commons_pc[A, B] > 0:
                        futures.append(executor.submit(calculate_significance, A, B, commons_pc, pcs, number_of_pc, logT, max_sig))

            for future in as_completed(futures):
                A, B, sig = future.result()
                if sig > thres:
                    S[A, B] = sig
                    S[B, A] = sig

        # 定期清理内存
        gc.collect()

    edge_count = np.count_nonzero(S > 0)
    min_sig = np.min(S[S > 0]) if edge_count > 0 else float('inf')
    max_sig = np.max(S[S > 0]) if edge_count > 0 else float('-inf')

    if edge_count > 0:
        print(f"Hypergeometric contig-similarity network:\n {contigs:10} contigs,\n {edge_count:10} edges (min: {min_sig:.2f}, max: {max_sig:.2f}, threshold was {thres})")
    else:
        raise ValueError("No edge in the similarity network!") 

    return S

# 加载数据并创建网络
matrix_path = f"{args.outpath}/out/matrix.npy"
singletons_path = f"{args.outpath}/out/singletons.npy"
thres = 1
max_sig = 300
S = create_network(matrix_path, singletons_path, thres, max_sig)

with open(f"{args.outpath}/out/output_network.npz", "wb") as f:
    np.savez(f, S=S)
