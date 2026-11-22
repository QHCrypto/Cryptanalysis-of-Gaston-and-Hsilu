"""Check the exact 5-round kernel repair and its small polynomial certificate.

Consumes complete sums from certificate.cpp and kernel_probe.cpp. No partial
cube sum is used to claim a coefficient or rank. Hex vectors pack word 0 low.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'src'))
from gaston_6r_ppm_transfer_audit import apply_columns, rotate_output
from structure import nullspace_from_columns, rank


def z_to_t(poly, width=64):
    """Substitute z=t+1 using Lucas' parity rule for binomial coefficients."""
    return sum((sum((poly >> j) & 1 for j in range(width) if j & i == i) & 1) << i
               for i in range(width))


def polynomial_product(a, b):
    product = 0
    while a:
        bit = a & -a
        product ^= b << (bit.bit_length() - 1)
        a ^= bit
    return product


def valuation(poly):
    return (poly & -poly).bit_length() - 1 if poly else None


def main():
    first = json.loads((HERE / 'certificate_r5.json').read_text())
    second = json.loads((HERE / 'kernel_probe_r5.json').read_text())
    assert first['coverage_complete'] and second['coverage_complete']
    assert first['rounds'] == second['rounds'] == 5
    assert first['cube_dimension'] == second['cube_dimension'] == 31
    assert first['covered_points'] == second['count'] == 1 << 31
    assert second['begin'] == 0
    representatives = [int(v, 16) for v in first['representative_columns_hex']]
    columns = [rotate_output(v, j, 64) for v in representatives for j in range(64)]
    kernel = nullspace_from_columns(columns)
    assert all(apply_columns(columns, v) == 0 for v in kernel)
    keys = [int(a, 16) | (int(b, 16) << 64) for a, b in second['keys_hex']]
    assert keys == [0] + kernel
    sums = [int(a, 16) | (int(b, 16) << 64) for a, b in second['sums_hex']]
    images = [v ^ sums[0] for v in sums[1:]]
    independent_rows, bit_positions = [], []
    for j in range(128):
        row = sum(((v >> j) & 1) << a for a, v in enumerate(images))
        if rank(independent_rows + [row]) > len(independent_rows):
            independent_rows.append(row)
            bit_positions.append(j)
        if len(independent_rows) == len(kernel):
            break

    # A circulant block is multiplication by p(z) modulo z^64+1=t^64.
    # A determinant with valuation d<64 certifies rank exactly 128-d.
    mask = (1 << 64) - 1
    polys = [[z_to_t((v >> (64 * w)) & mask) for v in representatives]
             for w in range(2)]
    determinant = (polynomial_product(polys[0][0], polys[1][1]) ^
                   polynomial_product(polys[0][1], polys[1][0])) & mask
    determinant_order = valuation(determinant)
    assert determinant_order is not None and determinant_order < 64
    assert rank(columns) == 128 - determinant_order
    coefficient_matrices = [[[p >> i & 1 for p in row] for row in polys]
                            for i in range(3)]
    # The first three coefficient matrices already determine the t^3 term:
    # M0=0, so no unknown t^3 coefficient can contribute at order three.
    truncated = [[p & 7 for p in row] for row in polys]
    small_det = (polynomial_product(truncated[0][0], truncated[1][1]) ^
                 polynomial_product(truncated[0][1], truncated[1][0])) & 15
    assert small_det == 8

    result = {
        'rounds': 5, 'word_bits': 64,
        'first_cube': 'N1=0; N0 repeats an even-parity 32-bit word twice',
        'first_cube_dimension': 31, 'first_cube_rank': rank(columns),
        'kernel_basis_hex': [f'{v:032x}' for v in kernel],
        'second_cube': second['cube'], 'second_cube_dimension': 31,
        'kernel_images_hex': [f'{v:032x}' for v in images],
        'rank_gain': rank(images), 'joint_rank': rank(columns) + rank(images),
        'repair_minor': {
            'packed_output_bit_positions': bit_positions,
            'rows': [[row >> a & 1 for a in range(len(kernel))]
                     for row in independent_rows],
            'rank': rank(independent_rows),
        },
        'sources': ['certificate_r5.json', 'kernel_probe_r5.json'],
        'complete_enumeration_points_per_cube': 1 << 31,
        'polynomial_certificate': {
            'coefficient_matrices_M0_M1_M2': coefficient_matrices,
            'entry_valuations': [[valuation(p) for p in row] for row in polys],
            'determinant_valuation': determinant_order,
            'smith_exponents': [1, 2],
            'rank_from_determinant': 128 - determinant_order,
            'proof': 'Over F2[t]/t^64, Smith exponents s1,s2 give rank 128-s1-s2. If the determinant has valuation d<64, neither exponent is saturated and s1+s2=d. Here the 12 bits in M0,M1,M2 give det=t^3+O(t^4), hence rank125.',
        },
        'scope': 'Exact actual-word-size five-round result; no six-round coefficient is supplied.',
    }
    path = HERE / 'kernel_repair_r5.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
