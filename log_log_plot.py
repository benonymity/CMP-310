from time import time
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def selection_sort(arr):
  n=len(arr)
  # make a local copy to modify and sort:
  result = list(arr)

  # compute result[i]
  for i in range(n):
    # find the minimum element in all remaining
    min_index=i

    for j in range(i+1,n):
      if result[j] < result[min_index]:
                min_index=j

    # swap(result[min_index], result[j])
    result[i],result[min_index] = result[min_index],result[i]
    
  return result

def selection_sort_slow(x):
  n = len(x)
  # after we find the first n-1 minimums, the final is already in the
  # right place:
  for i in range(n-1):
    # get the index at which min(x[i:]) occurs
    # (for brevity, done more hackily than in the notes):
    _,ind_min = min([ (x[j],j) for j in range(i,n) ])
    # swap x[i] and x[ind_min]
    x[i],x[ind_min] = x[i],x[ind_min]

def merge_sort(x:list):
  n = len(x)
  if n <= 1:
    return x
  L = merge_sort(x[:n//2])
  R = merge_sort(x[n//2:])
  # merge:
  res = []
  i_L = i_R = 0
  while i_L < len(L) and i_R < len(R):
    v_L = L[i_L]
    v_R = R[i_R]
    if v_L < v_R:
      res.append(v_L)
      i_L += 1
    else:
      res.append(v_R)
      i_R += 1

  res.extend(L[i_L:])
  res.extend(R[i_R:])
  return res

# the log problem sizes to be timed:
log_n_s = np.arange(3, 12)
# the problem sizes to be timed:
n_s = 2**log_n_s

num_reps = 2**6
# benchmark each algorithm on many problems on each of our problem
# sizes. for each measurement, take the average over num_reps
# measurements to reduce noise:

# 1) what are the \Theta runtimes of each algorithm timed?
# sum of a list \in Theta(???)
# merge_sort \in Theta(???)
# selection_sort \in Theta(???)
# selection_sort_slow \in Theta(???)

# 2) denote x=log(n) and y=log(r(n)) using the \Theta(.) runtime for
# each algorithm. for each algorithm, what is y as a function of x?

for alg in [sum, merge_sort, selection_sort, selection_sort_slow]:
  print(alg)
  times_for_alg = []
  for n in n_s:
    # take several empirical benchmarks and average:
    avg_time_at_n = 0.
    # NOTE: if you don't have the tqdm package installed, you can
    # replace the j for loop with the following line:
    #for j in range(num_reps):
    for j in tqdm(range(num_reps)):

      x = np.random.uniform(0., 1., n)
      t1 = time()
      alg(x)
      t2 = time()
      elapsed = t2 - t1
      avg_time_at_n += elapsed
    avg_time_at_n /= num_reps

    times_for_alg.append(avg_time_at_n)

  # plot the runtime series for this algorithm:
  plt.plot(n_s, times_for_alg, label=f'{alg}')

# display our runtimes with a legend:
plt.legend()
# 3) use the log-log plot to try to discuss the relative runtimes of
# these algorithms.

# a) how does a constant speedup manifest?

# b) how does an asymptotically better \Theta(.) manifest (i.e., a
# speedup that is more than a constant)?

# c) if one algorithm has a runtime \in \Theta(n^k) and another has a
# runtime \in \Theta(n^k \log(n)), what happens to the gap between
# them in the log-log plot as n grows? and what happens to the rate at
# which the gap between them on the log-log plot changes?

# use log axes for both x and y:
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$n$')
plt.ylabel('$r(n)$ in seconds')
plt.show()
