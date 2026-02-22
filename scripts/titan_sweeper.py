#!/usr/bin/env python3
"""
titan_sweeper.py — Quadruplet search for Q(n) = n^47 - (n-1)^47.

Scans a range [n_start, n_end] for 4-constellations: consecutive integers
n, n+1, n+2, n+3 such that all four Q(m) are probable primes.

Usage:
    python titan_sweeper.py <n_start> <n_end> [--output FILE]

Requirements:
    pip install gmpy2

Performance (3.5 GHz, 500-digit numbers):
    ~0.4 ms per candidate n (with pre-sieve)
    ~12 ms per Miller-Rabin test (25 rounds)

Author: Ruqing Chen, GUT Geoservice Inc.
Date:   February 2026
"""

import sys
import time
import argparse
import gmpy2


def Q(n):
    """Compute Q(n) = n^47 - (n-1)^47 using GMP arbitrary precision."""
    n = gmpy2.mpz(n)
    return n**47 - (n - 1)**47


def is_prime(x):
    """Probable prime test: 25-round Miller-Rabin via GMP."""
    return gmpy2.is_prime(x, 25)


def build_sieve_table(q):
    """
    For sieving prime q, compute the set of residues r in [0, q)
    such that Q(r) ≡ 0 (mod q).
    """
    bad = set()
    for r in range(q):
        if Q(r) % q == 0:
            bad.add(r)
    return bad


def sweep(n_start, n_end, output_file=None):
    """
    Scan [n_start, n_end] for prime quadruplets of Q(n).

    Strategy:
    1. Pre-sieve: for small primes q, precompute residues where Q(n) ≡ 0 mod q.
       Skip any n where Q(n), Q(n+1), Q(n+2), or Q(n+3) has a small factor.
    2. Primality: test surviving candidates with 25-round Miller-Rabin.
    3. Chain detection: maintain a sliding window of 4 consecutive prime flags.
    """
    # Build sieve tables for primes up to 1000
    sieve_primes = [p for p in range(2, 1000) if gmpy2.is_prime(p)]
    sieve_tables = {}
    print(f"Building sieve tables for {len(sieve_primes)} primes...")
    for q in sieve_primes:
        sieve_tables[q] = build_sieve_table(q)

    # Sliding window: track primality of Q(n), Q(n+1), Q(n+2), Q(n+3)
    quads_found = []
    quints_found = []
    tested = 0
    primes_found = 0
    t0 = time.time()

    # Pre-fill the window with n_start-3 .. n_start-1
    window = [False, False, False]  # Q(n_start-3), Q(n_start-2), Q(n_start-1)

    for n in range(n_start, n_end + 1):
        # Pre-sieve: check if Q(n) is divisible by any small prime
        sieved_out = False
        for q, bad_residues in sieve_tables.items():
            if (n % q) in bad_residues:
                sieved_out = True
                break

        if sieved_out:
            is_p = False
        else:
            val = Q(n)
            is_p = is_prime(val)
            tested += 1
            if is_p:
                primes_found += 1

        # Shift window
        window.append(is_p)

        # Check for quadruplet: window[-4:] = [True, True, True, True]
        if len(window) >= 4 and all(window[-4:]):
            start_n = n - 3
            digits = int(46 * gmpy2.log10(gmpy2.mpz(start_n)) + 1.67)
            quads_found.append((start_n, digits))
            print(f"  ★ QUADRUPLET at n = {start_n:,}  ({digits} digits)")

            # Check for quintuplet
            if len(window) >= 5 and window[-5]:
                quints_found.append(start_n - 1)
                print(f"  ★★ QUINTUPLET at n = {start_n - 1:,}")

        # Progress
        if (n - n_start) % 10_000_000 == 0 and n > n_start:
            elapsed = time.time() - t0
            pct = (n - n_start) / (n_end - n_start) * 100
            rate = (n - n_start) / elapsed
            print(f"  Progress: {pct:.1f}% | n = {n:,} | "
                  f"{rate:.0f} n/sec | {len(quads_found)} quads")

    elapsed = time.time() - t0

    # Report
    print(f"\n{'='*60}")
    print(f"Sweep complete: [{n_start:,}, {n_end:,}]")
    print(f"  Time: {elapsed:.1f} s ({elapsed/3600:.2f} h)")
    print(f"  Tested (post-sieve): {tested:,}")
    print(f"  Primes found: {primes_found:,}")
    print(f"  Quadruplets: {len(quads_found)}")
    print(f"  Quintuplets: {len(quints_found)}")
    print(f"{'='*60}")

    # Save results
    if output_file:
        with open(output_file, 'w') as f:
            f.write("starting_n,approx_digits\n")
            for start_n, digits in quads_found:
                f.write(f"{start_n},{digits}\n")
        print(f"Results saved to {output_file}")

    return quads_found, quints_found


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Quadruplet search for Q(n) = n^47 - (n-1)^47')
    parser.add_argument('n_start', type=int, help='Start of range')
    parser.add_argument('n_end', type=int, help='End of range')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output CSV file')
    args = parser.parse_args()

    sweep(args.n_start, args.n_end, args.output)
