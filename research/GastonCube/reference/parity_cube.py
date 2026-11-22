"""Full Boolean reference and rotation reconstruction for the even-parity cube."""
from gaston import IV, permute, rol


def nonce_word(index, rounds, width=64):
    period = 1 << rounds
    if period > width or width % period:
        raise ValueError('period 2^rounds must divide the word width')
    if index < 0 or index >= (1 << (period-1)):
        raise ValueError('cube assignment outside its range')
    block = (index << 1) | (index.bit_count() & 1)
    return sum(block << shift for shift in range(0, width, period))


def cube_sum(key, rounds, width=64, iv=IV, start=None):
    """Complete cube sum; exponential in 2^rounds-1. Intended for small checks."""
    result = 0
    mask = (1 << width)-1
    for index in range(1 << ((1 << rounds)-1)):
        state = [iv, key & mask, (key >> width) & mask,
                 nonce_word(index, rounds, width), 0]
        out = permute(state, rounds, width, start)
        result ^= out[0] | (out[1] << width)
    return result


def rotate_output(value, shift, width=64):
    mask = (1 << width)-1
    return rol(value & mask, shift, width) | (rol(value >> width, shift, width) << width)


def matrix_from_columns(columns, width=64):
    """Two packed representative columns, key bits 0 and width."""
    rotated = [rotate_output(columns[a], j, width) for a in range(2) for j in range(width)]
    return [[(col >> i) & 1 for col in rotated] for i in range(2*width)]


def matrix_from_rows(rows, width=64):
    """Rows for output bits 0 and width; missing entries stay None."""
    return [[rows[b][a*width+(j-t) % width]
             for a in range(2) for j in range(width)]
            for b in range(2) for t in range(width)]


def gf2_rank(matrix):
    pivots = {}
    for row in matrix:
        if any(x not in (0, 1) for x in row):
            raise ValueError('rank requires a fully certified binary matrix')
        value = sum(x << j for j, x in enumerate(row))
        while value:
            p = value.bit_length()-1
            if p not in pivots:
                pivots[p] = value
                break
            value ^= pivots[p]
    return len(pivots)
