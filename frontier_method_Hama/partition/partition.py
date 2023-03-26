def part_int_sub(n, k, a):
    if n == 0:   print(*a)
    elif n == 1: print(*(a + [1]))
    elif k == 1: print(*(a + [1] * n))
    # elif l == 10: 
    else:
        if n >= k:
            part_int_sub(n - k, k, a + [k])
        part_int_sub(n, k - 1, a)

def partition_of_int(n): part_int_sub(n, n, [])

part_int_sub(3, 3, [])


