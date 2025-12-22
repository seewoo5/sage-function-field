from functools import lru_cache


@lru_cache(maxsize=None)
def fibo_poly(n, q=None):
    """Fibonacci polynomial F_n over ZZ or GF(q)."""
    if q is None:
        R.<t> = ZZ['t']
    else:
        R.<t> = GF(q)['t']
    if n == 0:
        return R(0)
    elif n == 1:
        return R(1)
    else:
        return t * fibo_poly(n - 1, q) + fibo_poly(n - 2, q)


@lru_cache(maxsize=None)
def lucas_poly(n, q=None):
    """Lucas polynomial L_n over ZZ or GF(q)."""
    if q is None:
        R.<t> = ZZ['t']
    else:
        R.<t> = GF(q)['t']
    if n == 0:
        return R(2)
    elif n == 1:
        return t
    else:
        return t * lucas_poly(n - 1, q) + lucas_poly(n - 2, q)


@lru_cache(maxsize=None)
def W(n, f, g, q=None):
    """Generalized Lucas polynomial sequence W_n over ZZ or GF(q)."""
    if q is None:
        R.<t> = ZZ['t']
    else:
        R.<t> = GF(q)['t']
    if n == 0:
        return R(0)
    elif n == 1:
        return R(1)
    else:
        return f * W(n - 1, f, g, q) + g * W(n - 2, f, g, q)


@lru_cache(maxsize=None)
def w(n, f, g, q=None):
    """Generalized Lucas polynomial sequence w_n over ZZ or GF(q)."""
    if q is None:
        R.<t> = ZZ['t']
    else:
        R.<t> = GF(q)['t']
    if n == 0:
        return R(2)
    elif n == 1:
        return R(f)
    else:
        return f * w(n - 1, f, g, q) + g * w(n - 2, f, g, q)


@lru_cache(maxsize=None)
def U(n, q=None):
    """
    Chebyshev polynomial of the second kind U_n over ZZ or GF(q).
    U_n(T) = W_{n+1}(T) with f(T) = 2T and g(T) = -1.
    """
    if q is None:
        R.<t> = ZZ['t']
    else:
        R.<t> = GF(q)['t']
    return W(n + 1, 2 * t, R(-1), q)


@lru_cache(maxsize=None)
def J(n, q=None):
    """
    Jacobsthal polynomial J_n over ZZ or GF(q).
    J_n(T) = W_n(T) with f(T) = 1 and g(T) = 2T.
    """
    if q is None:
        R.<t> = ZZ['t']
    else:
        R.<t> = GF(q)['t']
    return W(n, R(1), 2 * t, q)


@lru_cache(maxsize=None)
def T(n, q=None):
    """
    Chebyshev polynomial of the first kind T_n over ZZ or GF(q).
    T_n(T) = w_n(T) / 2 with f(T) = 2T and g(T) = -1.
    """
    assert q is None or q % 2 == 1, "Only consider zero or odd characteristic for T_n."
    if q is None:
        R.<t> = ZZ['t']
        return w(n, 2 * t, R(-1), q) // 2
    else:
        R.<t> = GF(q)['t']
        return w(n, 2 * t, R(-1), q) * GF(q)(2)^(-1)

# For writing LaTeX tables. Written by ChatGPT-5.2.
def _poly_to_latex_with_var(poly, T_symbol="T"):
    s = poly._latex_()
    return s.replace("{t}", T_symbol).replace("t", T_symbol)


def _is_monomial_T(poly):
    """
    Check whether poly == t in GF(p)[t].
    """
    return poly.degree() == 1 and poly[0] == 0 and poly[1] == 1


def _latex_factorization_over_Fp(f, T_symbol="T", use_cdot=False):
    """
    Convert Sage factorization to LaTeX:
      T^{e} for pure powers of T,
      (factor)^{e} for all other factors.
    """
    fac = f.factor()
    parts = []

    u = f.lc()  # leading coefficient
    u_tex = _poly_to_latex_with_var(u, T_symbol=T_symbol)
    if u != 1:
        parts.append(u_tex)

    for g, e in fac:
        # Special case: pure power of T
        if _is_monomial_T(g):
            if e == 1:
                parts.append(T_symbol)
            else:
                parts.append(f"{T_symbol}^{{{e}}}")
            continue

        # General case
        gtex = _poly_to_latex_with_var(g, T_symbol=T_symbol)
        gtex = f"({gtex})"
        if e == 1:
            parts.append(gtex)
        else:
            parts.append(f"{gtex}^{{{e}}}")

    if not parts:
        return _poly_to_latex_with_var(f, T_symbol=T_symbol)

    joiner = " \\cdot " if use_cdot else ""
    return joiner.join(parts)


def print_latex_factor_table(
    N,
    p,
    poly_fn,
    poly_name="P",
    T_symbol="T",
    float_env="table",
    factor_col_width=None,   # e.g. "0.28\\textwidth" for each factor column
    fontsize_cmd=None,       # e.g. "\\scriptsize"
    caption=None,
    label=None,
    use_cdot=False,
):
    """
    Print a LaTeX table with columns:
      n | P_n(T) over F_p1 | P_n(T) over F_p2 | ...

    Parameters
    ----------
    N : int
        Print n = 1..N.
    p : int or list/tuple of int
        Prime(s). If int, treated as [p]. If list, each becomes a column.
    poly_fn : callable
        poly_fn(n, p) should return polynomial in GF(p)[t], e.g. lucas_poly(n,p).
    factor_col_width : str or None
        If provided, use p{<width>} for each factor column to allow wrapping.
        If None, uses l columns (no wrapping).
    """
    ps = [p] if isinstance(p, (int, Integer)) else list(p)

    # Column specification
    if factor_col_width is None:
        # no wrapping
        colspec = "c|" + "|".join(["l"] * len(ps))
    else:
        colspec = "c|" + "|".join([f"p{{{factor_col_width}}}"] * len(ps))

    # Header row
    header_cells = [f"$n$"] + [f"${poly_name}_n({T_symbol})$ over $\\mathbb{{F}}_{pp}$" for pp in ps]
    header = " & ".join(header_cells) + " \\\\"

    lines = []
    if float_env:
        lines.append(f"\\begin{{{float_env}}}")
        if fontsize_cmd:
            lines.append(fontsize_cmd)
        lines.append("\\centering")

    lines.append(f"\\begin{{tabular}}{{{colspec}}}")
    lines.append("\\toprule")
    lines.append(header)
    lines.append("\\midrule")

    for n in range(1, N + 1):
        row = [str(n)]
        for pp in ps:
            fnp = poly_fn(n, pp)
            tex = _latex_factorization_over_Fp(
                fnp,
                T_symbol=T_symbol,
                use_cdot=use_cdot,
            )
            row.append(f"${tex}$")
        lines.append(" & ".join(row) + " \\\\")

    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")

    if caption:
        lines.append(f"\\caption{{{caption}}}")
    if label:
        lines.append(f"\\label{{{label}}}")

    if float_env:
        lines.append(f"\\end{{{float_env}}}")

    print("\n".join(lines))
