#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The two decks' number formats, ported from their build_slides.js so
#           that a number printed in the workbook reads exactly as on the slide.
#           JavaScript's toFixed / toPrecision round the exact binary value half
#           up (0.125 -> "0.13"); Python's round() and format() round half to even
#           (0.125 -> "0.12"). Every format here therefore goes through Decimal on
#           the float's exact value with ROUND_HALF_UP.
#             pv     P: two significant digits from 0.001, else m.mm x 10^e (deck pv)
#             pe     P with two significant digits, for compact ranges (deck pe)
#             pspan  a P range on one exponent where possible (deck pspan)
#             fix    fixed decimals with a typographic minus (deck fix)
#             pct    100 x, fixed decimals, '%' (deck pct)
#             sig3   OR-style: 0 dp from 100, 1 dp from 10, else 2 dp (deck sig3)
#             integer  rounded, thousands separators (deck int)
#           Superscript exponents are returned as Unicode characters; the writer
#           turns them into superscript runs.
# Component: hla_records
# ---------------------------------------------------------------------------
from decimal import ROUND_HALF_UP, Decimal
from math import floor, isfinite, log10

SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
# V8 (node 22) evaluates 10 ** e for these exponents to a value one ulp off the exact power of ten; the decks'
# pv divides by it, so 9.995e-4 prints as 1.00 × 10⁻³ on the slide. Measured with node; elsewhere 10 ** e is exact.
V8_POW10 = {-29: 1.0000000000000001e-29, -24: 1.0000000000000001e-24, -23: 1.0000000000000001e-23,
            -21: 1.0000000000000001e-21, -20: 1.0000000000000001e-20, -17: 9.9999999999999992e-18,
            -5: 9.9999999999999991e-06, -4: 9.9999999999999991e-05}


def pow10(e):
    return V8_POW10.get(e, float(f'1e{e}'))
UNSUP = str.maketrans('⁻⁰¹²³⁴⁵⁶⁷⁸⁹', '-0123456789')


def to_fixed(x, dp):
    """JavaScript Number.prototype.toFixed: the exact binary value, rounded half up."""
    q = Decimal(1).scaleb(-dp) if dp > 0 else Decimal(1)
    d = Decimal(float(x)).quantize(q, rounding=ROUND_HALF_UP)
    s = f'{d:f}'
    if s.startswith('-') and Decimal(s) == 0:
        s = s[1:]
    return s


def to_precision(x, sd):
    """JavaScript Number.prototype.toPrecision for 1e-6 <= |x| < 1e21 (fixed notation)."""
    x = float(x)
    if x == 0:
        return to_fixed(0, sd - 1)
    e = floor(log10(abs(x)))
    d = Decimal(x).scaleb(-e).quantize(Decimal(1).scaleb(-(sd - 1)), rounding=ROUND_HALF_UP)
    if abs(d) >= 10:                     # 9.995 -> 10.0: one more order of magnitude
        e += 1
        d = Decimal(x).scaleb(-e).quantize(Decimal(1).scaleb(-(sd - 1)), rounding=ROUND_HALF_UP)
    v = d.scaleb(e)
    dp = max(sd - 1 - e, 0)
    return f'{v.quantize(Decimal(1).scaleb(-dp)) if dp else v.quantize(Decimal(1)):f}'


def sup(e):
    return str(e).translate(SUP)


def pv(p, dp=2):
    """The decks' P: p.toPrecision(2) from 0.001, else m (dp decimals) × 10^e."""
    p = float(p)
    if p >= 0.001:
        return to_precision(p, 2)
    e = floor(log10(p))
    m = p / pow10(e)
    if float(to_fixed(m, dp)) >= 10:
        m /= 10
        e += 1
    return f'{to_fixed(m, dp)} × 10{sup(e)}'


def pe(p):
    p = float(p)
    if p >= 0.01:
        return to_precision(p, 2)
    e = floor(log10(p))
    return f'{to_fixed(p / pow10(e), 1)} × 10{sup(e)}'


def pspan(xs):
    lo, hi = min(xs), max(xs)
    e = floor(log10(hi))
    if lo == hi:
        return pv(lo)
    if lo >= 0.001:
        return f'{pv(lo)}–{pv(hi)}'
    if hi >= 0.001 or floor(log10(lo)) != e:
        return f'{pv(lo)} to {pv(hi)}'
    a, b = to_fixed(lo / pow10(e), 2), to_fixed(hi / pow10(e), 2)
    return f'{a} × 10{sup(e)}' if a == b else f'{a}–{b} × 10{sup(e)}'


def fix(x, dp):
    return to_fixed(x, dp).replace('-', '−')


def pct(x, dp=1):
    return f'{to_fixed(100 * float(x), dp)}%'


def sig3(x):
    x = float(x)
    return to_fixed(x, 0 if x >= 100 else 1 if x >= 10 else 2)


def integer(x):
    return f'{int(Decimal(float(x)).quantize(Decimal(1), rounding=ROUND_HALF_UP)):,}'


def parse_p(text):
    """A displayed P back to a number: '4.57 × 10⁻⁵' -> 4.57e-5, '0.0026' -> 0.0026."""
    t = str(text).strip().translate(UNSUP).replace('−', '-')
    if '×' in t:
        m, e = t.split('×')
        return float(m) * 10 ** int(e.strip()[2:])
    return float(t)


def selftest():
    """Rounding cases where JavaScript and Python differ, and the deck's own examples."""
    cases = [(to_fixed(0.125, 2), '0.13'), (to_fixed(2.125, 2), '2.13'), (to_fixed(1.005, 2), '1.00'),
             (to_precision(0.00999, 2), '0.010'), (to_precision(0.0026, 2), '0.0026'), (pv(4.56965e-5), '4.57 × 10⁻⁵'),
             (pv(9.995e-4), '1.00 × 10⁻³'), (pv(0.104736), '0.10'), (pv(1.0), '1.0'), (pe(6.94e-3), '6.9 × 10⁻³'),
             (pct(0.10247, 2), '10.25%'), (sig3(0.65589), '0.66'), (integer(2661.5), '2,662'), (fix(-0.004, 2), '0.00'),
             (pspan([9.77e-5, 1.01e-4]), '9.77 × 10⁻⁵ to 1.01 × 10⁻⁴'), (pspan([2.48e-4, 2.49e-4]), '2.48–2.49 × 10⁻⁴')]
    bad = [(got, want) for got, want in cases if got != want]
    if bad:
        raise SystemExit(f'deckfmt selftest failed: {bad}')
    for s in ('4.57 × 10⁻⁵', '0.0026', '1.00 × 10⁻³'):
        assert abs(parse_p(s) - float(s.translate(UNSUP).replace(' × 10', 'e').replace('−', '-'))) < 1e-15 \
            if '×' in s else parse_p(s) == float(s)
    return len(cases)


if __name__ == '__main__':
    print(f'deckfmt: {selftest()} rounding cases agree with the decks')
