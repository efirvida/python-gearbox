"""Geometry factors ``I`` (pitting) and ``J`` (bending) of AGMA 908-B89 as used by AGMA 2001-D04."""

from __future__ import annotations

from dataclasses import dataclass
from math import acos, asin, atan, cos, degrees, pi, radians, sin, sqrt, tan

from gearbox.mathutils import inverse_involute
from gearbox.models import Gear
from gearbox.transmission import Transmission


@dataclass(frozen=True)
class GeometryFactors:
    """Geometry factors of a gear pair."""

    pitting: float
    """``I``."""
    bending: tuple[float, float]
    """``(J_pinion, J_wheel)``."""


@dataclass(frozen=True)
class _Mesh:
    """Mesh quantities shared by the ``J`` factor of both gears (AGMA 908-B89 symbols)."""

    psi: float
    fi_n: float
    fi_nr: float
    psi_r: float
    mf: float
    mn: float
    c_psi: float
    ro2: float


def geometry_factors(transmission: Transmission) -> GeometryFactors:
    """Compute the pitting and bending geometry factors, AGMA 908-B89 / AGMA 2001-D04.

    Glossary of the local variables (AGMA 908-B89 notation):

    * ``n1, n2``: numbers of teeth; ``nc``: teeth of the generating cutter.
    * ``mn``: normal module, later reused as the load sharing ratio ``m_N``.
    * ``fi_n, fi, fi_r, fi_nr, fi_nl``: normal, transverse, operating transverse, operating
      normal and load-application pressure angles; ``psi, psi_b, psi_r``: reference, base and
      operating helix angles. All in radians.
    * ``f``: face width over module; ``mp, mf``: transverse and face contact ratios.
    * ``ro1, ro2``: outside radii over module; ``r1, r2``: pitch radii; ``rb1, rb2``: base radii;
      ``mg``: gear ratio; ``cr``: centre distance over module.
    * ``c1..c6``: lengths along the line of action; ``pb, pn``: transverse and normal base pitch.
    * ``rho_1, rho_2``: radii of curvature at the point of load; ``c_psi``: helical overlap
      factor; ``lmin``: minimum length of the contact lines.
    * ``rn, rnb, n``: virtual spur gear pitch radius, base radius and teeth.
    * ``hf, sf, rho_f``: tooth root height, thickness and fillet radius of the critical section;
      ``kf``: stress correction factor; ``ch``: helical factor.

    Args:
        transmission: The gear pair.

    Returns:
        The :class:`GeometryFactors`.
    """
    p, w = transmission.pinion, transmission.wheel
    n1 = p.number_of_teeth
    n2 = w.number_of_teeth

    mn = p.module
    fi_n = radians(p.pressure_angle)
    psi = radians(p.helix_angle)
    f = p.face_width / mn

    ro1 = (p.tip_diameter / 2) / mn
    ro2 = (w.tip_diameter / 2) / mn
    mg = transmission.gear_ratio
    r1 = n1 / (2 * cos(psi))
    r2 = r1 * mg
    cr = r1 + r2
    fi = atan(tan(fi_n) / cos(psi))
    rb1 = r1 * cos(fi)
    rb2 = rb1 * mg

    if n1 > 0:
        fi_r = acos((rb2 + rb1) / cr)
    else:
        fi_r = acos((rb2 - rb1) / cr)

    pb = 2 * pi * rb1 / n1
    pn = pi * cos(fi_n)
    psi_b = acos(pn / pb)

    c6 = cr * sin(fi_r)
    if n1 > 0:
        c1 = c6 - sqrt(ro2**2 - rb2**2)
    else:
        c1 = -1 * (c6 - sqrt(ro2**2 - rb2**2))

    c4 = c1 + pb
    c5 = sqrt(ro1**2.0 - rb1**2.0)
    c2 = c5 - pb

    z = c5 - c1
    mp = z / pb

    nr = mp % 1
    if not psi:
        mf = 0
        mn = 1
        lmin = f
    else:
        px = pi / sin(psi)
        mf = f / px
        if mf > 1:
            na = mf % 1
            if 1 - nr < na:
                lmin = (mp * f - (1.0 - na) * (1.0 - nr) * px) / cos(psi_b)
            else:
                lmin = (mp * f - na * nr * px) / cos(psi_b)
            mn = f / lmin
        else:
            mn = 1
    psi_r = atan(tan(psi_b) / cos(fi_r))
    fi_nr = asin(cos(psi_b) * sin(fi_r))

    # I factor
    if n1 > 0:
        d = 2 * cr / (mg + 1.0)
        rm1 = 0.5 * (ro1 + (cr - ro2))
    else:
        d = 2 * cr / (mg - 1.0)
        rm1 = 0.5 * (ro1 + (cr - ro2))

    if mf <= 1:
        rho_1 = c2
        rho_m1 = sqrt(rm1**2 - rb1**2)
        if n1 > 0:
            rho_2 = c6 - rho_1
            rho_m2 = c6 - rho_m1
        else:
            rho_2 = c6 + rho_1
            rho_m2 = c6 + rho_m1
        c_psi = sqrt(1 - mf * (1 - rho_m1 * rho_m2 * z / (rho_1 * rho_2 * pn)))
    else:
        rho_1 = sqrt(rm1**2.0 - rb1**2.0)
        if n1 > 0:
            rho_2 = c6 - rho_1
        else:
            rho_2 = c6 + rho_1
        c_psi = 1

    if n1 > 0:
        pitting = (cos(fi_r) * c_psi**2.0) / ((1.0 / rho_1 + 1.0 / rho_2) * d * mn)
    else:
        pitting = (cos(fi_r) * c_psi**2.0) / ((1.0 / rho_1 - 1.0 / rho_2) * d * mn)

    # KNOWN-ISSUE KI-15: ``mn`` (reused above as the load sharing ratio) and ``ro2`` are shared
    # by both J factors, so the wheel J factor is evaluated with the pinion outside radius.
    mesh = _Mesh(
        psi=psi,
        fi_n=fi_n,
        fi_nr=fi_nr,
        psi_r=psi_r,
        mf=mf,
        mn=mn,
        c_psi=c_psi,
        ro2=ro2,
    )
    j1 = _bending_geometry_factor(
        mesh, p, n1, mg, ro1, r1, r2, rb1, c4, p.profile_shift_coefficient, p.backlash
    )
    j2 = _bending_geometry_factor(
        mesh, w, n2, 1.0 / mg, ro2, r2, r1, rb2, c6 - c2, w.profile_shift_coefficient, w.backlash
    )
    return GeometryFactors(pitting=pitting, bending=(j1, j2))


def _bending_geometry_factor(
    mesh: _Mesh,
    gear: Gear,
    n1: int,
    mg: float,
    ro1: float,
    r1: float,
    r2: float,
    rb1: float,
    c4: float,
    x: float,
    delta_sn: float,
) -> float:
    """Bending strength geometry factor ``J`` of one gear (AGMA 908-B89, clause 8).

    Iterates on the load angle ``alpha_n`` until the Lewis parabola is tangent to the root
    fillet. Symbols are those of the standard (see :func:`geometry_factors`); additionally
    ``xg`` is the profile shift corrected by the backlash, ``sn`` the tooth thickness,
    ``rno``/``rnbo``: virtual cutter pitch and base radius, ``r_s_no``: radius of the cutter
    tip centre, ``mi_no``, ``sigma_n``, ``beta_n``, ``xi_nf``, ``eta_nf``: auxiliary angles
    and coordinates of the fillet point, and ``y``/``y_`` the residual of the iteration and
    its derivative.

    Args:
        mesh: Quantities shared by both gears.
        gear: The gear whose factor is computed (supplies the tool data).
        n1: Number of teeth of the gear.
        mg: Gear ratio seen from this gear (1 for the pinion).
        ro1: Outside radius over module.
        r1: Pitch radius over module.
        r2: Pitch radius over module of the mate.
        rb1: Base radius over module.
        c4: Length along the line of action to the lowest point of single tooth contact.
        x: Profile shift coefficient.
        delta_sn: Backlash allowance.

    Returns:
        The geometry factor ``J``.
    """
    psi, fi_n, fi_nr, psi_r = mesh.psi, mesh.fi_n, mesh.fi_nr, mesh.psi_r
    mf, mn, c_psi, ro2 = mesh.mf, mesh.mn, mesh.c_psi, mesh.ro2
    tool = gear.tool
    nc = tool.cutter_teeth
    hao = tool.addendum_coefficient
    xo = tool.profile_shift_coefficient
    rho_ao = tool.tool_edge_radius_coefficient
    delta_ao = tool.tool_tip_displacement

    if psi != 0:
        n = n1 / cos(psi) ** 3.0
        rn = n / 2
        rnb = rn * cos(fi_n)
        if mf <= 1:
            rn2 = rn * mg
            rnb2 = rnb * mg
            rna2 = rn2 + ro2 - r2
            cn6 = (rnb2 + rnb) * tan(fi_nr)
            cn1 = cn6 - sqrt(rna2**2 - rnb2**2)
            cn4 = cn1 + pi * cos(fi_n)
            tan_fi_nw = cn4 / rnb
        else:
            rna = rn + ro1 - r1
            tan_fi_nw = sqrt(((rna / rnb) ** 2.0) - 1.0)
    else:
        n = n1
        rn = r1
        rnb = rb1
        tan_fi_nw = c4 / rnb

    xg = x - delta_sn / (2 * tan(fi_n))
    sn = pi / 2 + 2 * xg * tan(fi_n)
    inv_fi_n = tan(fi_n) - fi_n
    fi_nl = tan_fi_nw - tan(fi_n) + fi_n - sn / n
    rnl = rnb / cos(fi_nl)
    no = nc / cos(psi) ** 3.0
    rno = no / 2
    rnbo = rno * cos(fi_n)
    r_s_no = rno + hao + xo - rho_ao
    fi_ns = acos(rnbo / r_s_no)
    inv_fi_ns = tan(fi_ns) - fi_ns
    sno = pi / 2 + 2 * xo * tan(fi_n)
    inv_fi_npo = tan(fi_n) - fi_n + sno / no
    lambda_ns_2 = inv_fi_npo - inv_fi_ns + ((delta_ao - rho_ao) / rnbo)
    inv_fi_n_2 = inv_fi_n + ((2 * (xg + xo) * tan(fi_n)) / (n + no))
    fi_n_2 = radians(inverse_involute(inv_fi_n_2))

    rn_2 = rn * cos(fi_n) / cos(fi_n_2)
    rno_2 = rno * cos(fi_n) / cos(fi_n_2)
    alpha_n = pi / 4
    while True:
        mi_no = acos(rno_2 * cos(alpha_n) / r_s_no) - alpha_n
        ks = rno_2 * sin(alpha_n) - r_s_no * sin(alpha_n + mi_no)
        kf = ks - rho_ao
        sigma_n = no / n * (mi_no - lambda_ns_2 + pi / no)
        beta_n = alpha_n - sigma_n
        xi_nf = rn_2 * sin(sigma_n) + kf * cos(beta_n)
        eta_nf = rn_2 * cos(sigma_n) + kf * sin(beta_n)
        hf = rnl - eta_nf
        y = 2 * hf * tan(beta_n) - xi_nf
        y_ = (
            (2 * hf / cos(beta_n) ** 2)
            - kf * sin(beta_n)
            + no
            / n
            * ((rno_2 * sin(alpha_n) / (r_s_no * sin(alpha_n + mi_no))) - 1)
            * (2 * xi_nf * tan(beta_n) - eta_nf - 2 * hf / (cos(beta_n) ** 2))
            - rno_2
            * (cos(alpha_n) - sin(alpha_n) / tan(alpha_n + mi_no))
            * ((1 + sin(beta_n) ** 2.0) / cos(beta_n))
        )
        alpha_n1 = alpha_n - y / y_

        # KNOWN-ISSUE KI-15: the convergence test is on the signed residual (``y <= 1e-6``), not
        # its absolute value, so the loop may stop on a negative, unconverged residual.
        if y <= 1e-6:
            break
        alpha_n = alpha_n1

    rho_f = rho_ao + ((rno_2 - r_s_no) ** 2) / ((rn_2 * rno_2 / (rn_2 + rno_2)) - (rno_2 - r_s_no))
    omega = degrees(atan(tan(psi) * sin(fi_n)))
    sf = 2.0 * xi_nf
    h = 0.331 - 0.436 * fi_n
    l = 0.324 - 0.492 * fi_n  # noqa: E741
    m = 0.261 - 0.545 * fi_n
    kf = h + ((sf / rho_f) ** l) * ((sf / hf) ** m)

    if mf > 1:
        ch = 1.0 / (1.0 - sqrt((omega / 100.0) * (1.0 - omega / 100.0)))
        k_psi = cos(psi_r) * cos(psi)
    else:
        ch = 1.0
        k_psi = 1.0

    y = k_psi / ((cos(fi_nl) / cos(fi_nr)) * ((6.0 * hf / ((sf**2.0) * ch)) - tan(fi_nl) / sf))

    # KNOWN-ISSUE KI-15: for the pinion of the "spur_sr0"/"spur_sr20" scenarios J exceeds 1
    # (1.045); AGMA 908-B89 J factors of external spur gears are normally below 1.
    return y * c_psi / (kf * mn)
