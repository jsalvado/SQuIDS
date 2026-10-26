#!/usr/bin/env python3
"""Generate include/SQuIDS/detail/rk_tableaus.h.

Coefficients are parsed from the published sources listed in SOURCES, kept as
exact rationals (Fraction), checked against every Runge-Kutta order condition
(all rooted trees up to the claimed order), and only then rounded to double.

    python3 tools/gen_rk_tableaus.py <dir with downloaded sources> > include/SQuIDS/detail/rk_tableaus.h

The source files are not stored in the repository; SOURCES records where each
one came from.
"""
import re
import sys
from fractions import Fraction
from functools import reduce

SOURCES = {
    "verner_tableaus.jl": "https://github.com/SciML/OrdinaryDiffEq.jl lib/OrdinaryDiffEqVerner/src/verner_tableaus.jl",
    "tsit_tableaus.jl": "https://github.com/SciML/OrdinaryDiffEq.jl lib/OrdinaryDiffEqTsit5/src/tsit_tableaus.jl",
    "rk.py": "https://github.com/scipy/scipy scipy/integrate/_ivp/rk.py",
    "dop853_coefficients.py": "https://github.com/scipy/scipy scipy/integrate/_ivp/dop853_coefficients.py",
}

# --------------------------------------------------------------------------
# Exact arithmetic in Q(sqrt 6), needed for Verner's 9(8) pair


class Q6:
    """a + b*sqrt(6) with rational a, b."""

    def __init__(self, a, b=0):
        self.a, self.b = Fraction(a), Fraction(b)

    @staticmethod
    def lift(x):
        return x if isinstance(x, Q6) else Q6(x)

    def __add__(self, o):
        o = Q6.lift(o)
        return Q6(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return Q6(-self.a, -self.b)

    def __sub__(self, o):
        return self + (-Q6.lift(o))

    def __rsub__(self, o):
        return Q6.lift(o) - self

    def __mul__(self, o):
        o = Q6.lift(o)
        return Q6(self.a * o.a + 6 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = Q6.lift(o)
        n = o.a * o.a - 6 * o.b * o.b
        return self * Q6(o.a / n, -o.b / n)

    def __rtruediv__(self, o):
        return Q6.lift(o) / self

    def __eq__(self, o):
        o = Q6.lift(o)
        return self.a == o.a and self.b == o.b

    def __ne__(self, o):
        return not self == o

    def __bool__(self):
        return self.a != 0 or self.b != 0

    def __float__(self):
        from decimal import Decimal, getcontext
        getcontext().prec = 80
        d = lambda f: Decimal(f.numerator) / Decimal(f.denominator)
        return float(str(d(self.a) + d(self.b) * Decimal(6).sqrt()))

    def __abs__(self):
        return abs(float(self))


def simplify(x):
    return x.a if isinstance(x, Q6) and x.b == 0 else x


# --------------------------------------------------------------------------
# Expression parsing


def to_fraction(expr):
    """Evaluate a Julia/Python coefficient expression exactly."""
    expr = " ".join(expr.split())
    lits = {}

    def stash(m):
        key = "L%d" % len(lits)
        lits[key] = Fraction(m.group(1))
        return key

    # high-precision decimal literals in their various spellings
    expr = re.sub(r'big"([^"]+)"', stash, expr)
    expr = re.sub(r'BigFloat\(\s*"([^"]+)"\s*\)', stash, expr)
    expr = re.sub(r'parse\(\s*BigFloat\s*,\s*"([^"]+)"\s*\)', stash, expr)
    expr = expr.replace("BigInt", "")
    # sqrt(6), written as 6^(1 // 2) in the Verner 9(8) source
    expr = re.sub(r"6\s*\^\s*\(\s*1\s*//\s*2\s*\)", "S6", expr)
    lits["S6"] = Q6(0, 1)
    # plain decimal literals, kept exact as written (only in non-rational sources)
    expr = re.sub(r"(?<![\w.])(\d+\.\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?|\d+[eE][-+]?\d+)",
                  stash, expr)
    expr = expr.replace("//", "/")
    expr = re.sub(r"(?<![\w.])(\d+)(?![\w.])", r"Fraction(\1)", expr)
    val = eval(expr, {"Fraction": Fraction, "__builtins__": {}}, lits)
    return simplify(val) if isinstance(val, Q6) else Fraction(val)


def balanced(text, start):
    """Return the contents of the parenthesis opened at text[start]."""
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                return text[start + 1:i]
    raise ValueError("unbalanced")


def julia_constructor(path, signature):
    """Parse `name = convert(T, expr)` assignments in a generic constructor."""
    src = open(path).read()
    i = src.index(signature)
    j = src.index("\nend", i)
    body = src[i:j]
    out = {}
    for m in re.finditer(r"(\w+)\s*=\s*convert\(", body):
        arg = balanced(body, m.end() - 1)
        _, expr = arg.split(",", 1)
        out[m.group(1)] = to_fraction(expr)
    return out


# --------------------------------------------------------------------------
# Tableau container


class Tableau:
    def __init__(self, name, order, err_order, A, b, e, fsal, source, note=""):
        self.name, self.order, self.err_order = name, order, err_order
        self.s = len(b)
        self.A = A          # s x s, strictly lower triangular
        self.b = b          # propagated solution weights
        self.e = e          # error weights: err = h * sum e_i k_i
        self.fsal = fsal    # last stage is f(t+h, y_new)
        self.c = [sum(row) for row in A]
        self.source, self.note = source, note


def square(s):
    return [[Fraction(0)] * s for _ in range(s)]


def from_names(coef, s, row_col, b_idx, e_idx, b_name="b", e_name="btilde", fsal_row=None):
    """Assemble a tableau from named coefficients.

    row_col(name) -> (row, col) 1-based for an `a...` name, or None.
    fsal_row: the A row (1-based) that holds the solution weights.
    """
    A = square(s)
    for k, v in coef.items():
        rc = row_col(k)
        if rc:
            A[rc[0] - 1][rc[1] - 1] = v
    if fsal_row:
        b = A[fsal_row - 1][:]
    else:
        b = [Fraction(0)] * s
        for i in b_idx:
            b[i - 1] = coef["%s%d" % (b_name, i)]
    e = [Fraction(0)] * s
    for i in e_idx:
        e[i - 1] = coef["%s%d" % (e_name, i)]
    return A, b, e


def rc_short(k):
    """a{i}{j} with single-digit indices (Vern6, Tsit5)."""
    m = re.fullmatch(r"a(\d)(\d)", k)
    return (int(m.group(1)), int(m.group(2))) if m else None


def rc_vern7(k):
    """a + 2-digit row + 1-digit column."""
    m = re.fullmatch(r"a(\d\d)(\d)", k)
    return (int(m.group(1)), int(m.group(2))) if m else None


def rc_wide(k):
    """a + 2-digit row + 2-digit column (Vern8, Vern9)."""
    m = re.fullmatch(r"a(\d\d)(\d\d)", k)
    return (int(m.group(1)), int(m.group(2))) if m else None


def rc_rkv76(k):
    """a{row}{col} with rows up to 10 and single-digit columns."""
    m = re.fullmatch(r"a(10|\d)(\d)", k)
    return (int(m.group(1)), int(m.group(2))) if m else None


# --------------------------------------------------------------------------
# Methods


def build(srcdir):
    out = []
    v = srcdir + "/verner_tableaus.jl"
    t5 = srcdir + "/tsit_tableaus.jl"

    # Fehlberg 4(5), exactly as GSL's rkf45.c (propagates the 5th-order solution)
    s = 6
    A = square(s)
    F = Fraction
    A[1][0] = F(1, 4)
    A[2][:2] = [F(3, 32), F(9, 32)]
    A[3][:3] = [F(1932, 2197), F(-7200, 2197), F(7296, 2197)]
    A[4][:4] = [F(439, 216), F(-8), F(3680, 513), F(-845, 4104)]
    A[5][:5] = [F(-8, 27), F(2), F(-3544, 2565), F(1859, 4104), F(-11, 40)]
    b5 = [F(16, 135), 0, F(6656, 12825), F(28561, 56430), F(-9, 50), F(2, 55)]
    b4 = [F(25, 216), 0, F(1408, 2565), F(2197, 4104), F(-1, 5), 0]
    out.append(Tableau("rkf45", 5, 4, A, [F(x) for x in b5], [F(x) - F(y) for x, y in zip(b5, b4)],
                       False, "Fehlberg (1969); identical to GSL ode-initval2/rkf45.c"))

    # Bogacki-Shampine 3(2), FSAL; and Dormand-Prince 5(4), FSAL (scipy RK23/RK45)
    out.append(scipy_rk(srcdir + "/rk.py", "RK23", "bs3", 3, 2))
    out.append(scipy_rk(srcdir + "/rk.py", "RK45", "dp5", 5, 4))

    # Tsitouras 5(4), FSAL
    c = julia_constructor(t5, "Tsit5ConstantCacheActual(::Type{T}, ::Type{T2}) where {T, T2}")
    A, b, e = from_names(c, 7, rc_short, None, range(1, 8), fsal_row=7)
    out.append(Tableau("tsit5", 5, 4, A, b, e, True, "Tsitouras (2011), " + SOURCES["tsit_tableaus.jl"]))

    # Verner "efficient" 6(5), FSAL
    c = julia_constructor(v, "Vern6Tableau(::Type{T}, ::Type{T2}) where {T, T2}")
    A, b, e = from_names(c, 9, rc_short, None, [1, 4, 5, 6, 7, 8, 9], fsal_row=9)
    out.append(Tableau("vern6", 6, 5, A, b, e, True, "Verner (2010), " + SOURCES["verner_tableaus.jl"]))

    # Verner "efficient" 7(6)
    c = julia_constructor(v, "Vern7Tableau(::Type{T}, ::Type{T2}) where {T, T2}")
    A, b, e = from_names(c, 10, rc_vern7, [1, 4, 5, 6, 7, 8, 9], [1, 4, 5, 6, 7, 8, 9, 10])
    out.append(Tableau("vern7", 7, 6, A, b, e, False, "Verner (2010), " + SOURCES["verner_tableaus.jl"]))

    # Verner "efficient" 8(7)
    c = julia_constructor(v, "Vern8Tableau(::Type{T}, ::Type{T2}) where {T, T2}")
    A, b, e = from_names(c, 13, rc_wide, [1, 6, 7, 8, 9, 10, 11, 12], [1, 6, 7, 8, 9, 10, 11, 12, 13])
    out.append(Tableau("vern8", 8, 7, A, b, e, False, "Verner (2010), " + SOURCES["verner_tableaus.jl"]))

    # Verner "efficient" 9(8)
    c = julia_constructor(v, "Vern9Tableau(::Type{T}, ::Type{T2}) where {T, T2}")
    A, b, e = from_names(c, 16, rc_wide, [1] + list(range(8, 16)), [1] + list(range(8, 17)))
    out.append(Tableau("vern9", 9, 8, A, b, e, False, "Verner (2010), " + SOURCES["verner_tableaus.jl"]))

    # Verner RKV76.IIa "most efficient" 7(6); error = b - bh
    c = julia_constructor(v, "RKV76IIaTableau(::Type{T}, ::Type{T2}) where {T, T2}")
    A = square(10)
    for k, val in c.items():
        rc = rc_rkv76(k)
        if rc:
            A[rc[0] - 1][rc[1] - 1] = val
    b = [c["b%d" % i] for i in range(1, 11)]
    bh = [c["bh%d" % i] for i in range(1, 11)]
    out.append(Tableau("rkv76", 7, 6, A, b, [x - y for x, y in zip(b, bh)], False,
                       "Verner RKV76.IIa, " + SOURCES["verner_tableaus.jl"]))

    out.append(dop853(srcdir + "/dop853_coefficients.py"))
    return out


def scipy_rk(path, cls, name, p, q):
    src = open(path).read()
    i = src.index("class %s(" % cls)
    j = src.index("\nclass ", i + 1)
    body = src[i:j]

    def arr(field):
        m = re.search(r"\n\s+%s = np\.array\(" % field, body)
        txt = balanced(body, m.end() - 1)
        rows = re.findall(r"\[([^\[\]]*)\]", txt)
        return [[to_fraction(x) for x in r.split(",") if x.strip()] for r in rows]

    A0 = arr("A")
    s0 = len(A0)
    B = arr("B")[0]
    E = arr("E")[0]
    s = s0 + 1  # scipy stores the FSAL stage implicitly (row s = B)
    A = square(s)
    for r in range(s0):
        for k in range(len(A0[r])):
            A[r][k] = A0[r][k]
    A[s0][:s0] = B
    b = B + [Fraction(0)]
    # scipy's E is (embedded - propagated); this file uses e = b - bhat
    return Tableau(name, p, q, A, b, [-x for x in E], True, "scipy %s, %s" % (cls, SOURCES["rk.py"]))


def dop853(path):
    src = open(path).read()
    n = 12
    A = square(16)
    for m in re.finditer(r"^A\[(\d+), (\d+)\] = (\S+)", src, re.M):
        A[int(m.group(1))][int(m.group(2))] = to_fraction(m.group(3))
    B = A[n][:n]
    E3 = B + [Fraction(0)]
    for m in re.finditer(r"^E3\[(\d+)\] -= (\S+)", src, re.M):
        E3[int(m.group(1))] -= to_fraction(m.group(2))
    E5 = [Fraction(0)] * (n + 1)
    for m in re.finditer(r"^E5\[(\d+)\] = (\S+)", src, re.M):
        E5[int(m.group(1))] = to_fraction(m.group(2))
    s = n + 1  # 12 stages + FSAL evaluation f(t+h, y_new)
    At = [row[:s] for row in A[:s]]
    t = Tableau("dop853", 8, 7, At, B + [Fraction(0)], E5, True,
                "Hairer DOP853 via " + SOURCES["dop853_coefficients.py"],
                note="error = |E5| * |E5|/sqrt(|E5|^2 + 0.01 |E3|^2)")
    t.e3 = E3
    return t


# --------------------------------------------------------------------------
# Order conditions


def rooted_trees(pmax):
    """All rooted trees up to order pmax as (children_indices, order, gamma)."""
    trees = [((), 1, 1)]
    by_order = {1: [0]}
    for n in range(2, pmax + 1):
        found = []

        def gen(rem, min_idx, acc):
            if rem == 0:
                found.append(tuple(acc))
                return
            for idx in range(min_idx, len(trees)):
                o = trees[idx][1]
                if o <= rem:
                    gen(rem - o, idx, acc + [idx])

        gen(n - 1, 0, [])
        by_order[n] = []
        for ch in found:
            g = n * reduce(lambda x, y: x * y, (trees[k][2] for k in ch), 1)
            by_order[n].append(len(trees))
            trees.append((ch, n, g))
    return trees


def check_order(t, weights, p, exact):
    trees = rooted_trees(p)
    s = t.s
    Ag = []  # A @ g(tree), memoized per tree
    worst = Fraction(0)
    for ch, o, gam in trees:
        v = [Fraction(1)] * s
        for k in ch:
            v = [x * y for x, y in zip(v, Ag[k])]
        Ag.append([sum(t.A[i][j] * v[j] for j in range(i) if t.A[i][j]) for i in range(s)])
        res = sum(w * x for w, x in zip(weights, v)) - Fraction(1, gam)
        worst = max(worst, abs(res))
        if exact and res != 0:
            return False, "tree of order %d fails exactly (residual %.3e)" % (o, float(res))
    if not exact and worst > 1e-14:
        return False, "max residual %.3e" % float(worst)
    return True, ("exact" if exact else "max residual %.1e" % float(worst))


def is_rational_source(t):
    return t.name in ("rkf45", "bs3", "dp5", "vern6", "vern7", "vern8", "vern9")


# --------------------------------------------------------------------------
# Output


def cdouble(x):
    return repr(float(x))


def emit(tabs, log):
    print("// Generated by tools/gen_rk_tableaus.py -- do not edit by hand.")
    print("// Every tableau below was checked against all Runge-Kutta order conditions")
    print("// (rooted trees) for both the propagated and the embedded solution.")
    print("#ifndef SQUIDS_DETAIL_RK_TABLEAUS_H")
    print("#define SQUIDS_DETAIL_RK_TABLEAUS_H")
    print("namespace squids{ namespace detail{ namespace rk{")
    print()
    for t in tabs:
        s = t.s
        print("// %s: order %d(%d), %d stages%s. %s" % (t.name, t.order, t.err_order, s,
              ", FSAL" if t.fsal else "", t.source))
        if t.note:
            print("// %s" % t.note)
        print("// order check: %s" % log[t.name])
        print("static const double %s_c[%d] = {%s};" % (t.name, s, ", ".join(cdouble(x) for x in t.c)))
        flat = [t.A[i][j] for i in range(s) for j in range(s)]
        print("static const double %s_A[%d] = {%s};" % (t.name, s * s, ", ".join(cdouble(x) for x in flat)))
        print("static const double %s_b[%d] = {%s};" % (t.name, s, ", ".join(cdouble(x) for x in t.b)))
        print("static const double %s_e[%d] = {%s};" % (t.name, s, ", ".join(cdouble(x) for x in t.e)))
        if hasattr(t, "e3"):
            print("static const double %s_e3[%d] = {%s};" % (t.name, s, ", ".join(cdouble(x) for x in t.e3)))
        print()
    print("}}} // namespace squids::detail::rk")
    print("#endif")


def main():
    srcdir = sys.argv[1]
    tabs = build(srcdir)
    log = {}
    for t in tabs:
        exact = is_rational_source(t)
        ok_b, msg_b = check_order(t, t.b, t.order, exact)
        bhat = [x - y for x, y in zip(t.b, t.e)]
        if t.name == "dop853":
            # two embedded solutions, of orders 5 and 3, combined into one estimate
            ok5, m5 = check_order(t, bhat, 5, exact)
            ok3, m3 = check_order(t, [x - y for x, y in zip(t.b, t.e3)], 3, exact)
            ok_e, msg_e = ok5 and ok3, "5th %s, 3rd %s" % (m5, m3)
        else:
            ok_e, msg_e = check_order(t, bhat, t.err_order, exact)
        # the error weights must not themselves be consistent (sum e = 0 is required)
        if sum(t.e) != 0 and exact:
            ok_e, msg_e = False, "sum of error weights is %s, expected 0" % sum(t.e)
        if t.fsal:
            assert t.A[-1][:-1] == t.b[:-1] and t.b[-1] == 0, t.name + ": FSAL row mismatch"
            off = abs(t.c[-1] - 1)
            assert (off == 0) if exact else (off < 1e-25), t.name + ": FSAL stage not at t+h"
        if not (ok_b and ok_e):
            sys.exit("%s FAILED: propagated %s; embedded %s" % (t.name, msg_b, msg_e))
        log[t.name] = "order %d %s; embedded order %d %s" % (t.order, msg_b, t.err_order, msg_e)
        sys.stderr.write("%-7s %s\n" % (t.name, log[t.name]))
    emit(tabs, log)


if __name__ == "__main__":
    main()
