from multiprocessing import Pool
pool = Pool()
result1 = pool.apply_async(solve1, [A])    # evaluate "solve1(A)" asynchronously
result2 = pool.apply_async(solve2, [B])    # evaluate "solve2(B)" asynchronously
answer1 = result1.get(timeout=10)
answer2 = result2.get(timeout=10)