"""Stereo70 and ETRS89 Geodetic Transformation Coordinator.

This module provides the TransRO class to coordinate 2D and 3D coordinate transformations
between Romania's official national cartographic projection (Stereo 70) and the
European Terrestrial Reference System 1989 (ETRS89).

Transformation Pipelines:
    - Forward (Stereo 70 -> ETRS89):
        1. 2D grid inverse shift correction (:class:`~pytransdatro.trans_grid.Grid2D`)
        2. 2D Helmert conformal transformation (:class:`~pytransdatro.trans_helmert2d.Helmert2D`)
        3. Oblique Stereographic projection to ellipsoidal coordinates (:class:`~pytransdatro.proj_stereo.StereoProj`)
        4. Optional 1D geoid height correction: :math:`h = H + N` (:class:`~pytransdatro.trans_grid.Grid1D`)

    - Inverse (ETRS89 -> Stereo 70):
        1. Oblique Stereographic conversion to projected coordinates (:class:`~pytransdatro.proj_stereo.StereoProj`)
        2. 2D Helmert inverse transformation (:class:`~pytransdatro.trans_helmert2d.Helmert2D`)
        3. 2D grid forward shift correction (:class:`~pytransdatro.trans_grid.Grid2D`)
        4. Optional 1D geoid height correction: :math:`H = h - N` (:class:`~pytransdatro.trans_grid.Grid1D`)
"""

import copy
import logging
import math
from typing import Any, Optional, Sequence, Tuple, Union, overload

from pytransdatro import exceptions
from pytransdatro.proj_stereo import StereoProj
from pytransdatro.trans_grid import Grid1D, Grid2D
from pytransdatro.trans_helmert2d import Helmert2D
from pytransdatro.types import Coord2D, Coord3D, ScalarOrSeq, Sequence2D, Sequence3D

logger = logging.getLogger(__name__)


def _validate_and_normalize_inputs(
    c1: ScalarOrSeq,
    c2: ScalarOrSeq,
    c3: Optional[ScalarOrSeq] = None,
) -> Tuple[bool, Sequence[float], Sequence[float], Optional[Sequence[float]]]:
    """Validates and normalizes scalar or sequence coordinate inputs.

    Enforces:
        1. Non-empty sequences when sequence inputs are provided.
        2. Matching dimensionality (scalar vs sequence) and equal sequence lengths.
        3. Finite numeric values (rejecting None, strings, NaN, Infinity).

    :param c1: First coordinate axis (Northing or Latitude)
    :type c1: float or Sequence[float]
    :param c2: Second coordinate axis (Easting or Longitude)
    :type c2: float or Sequence[float]
    :param c3: Optional third coordinate axis (Normal or Ellipsoidal height)
    :type c3: float or Sequence[float], optional
    :return: Tuple of ``(is_scalar, c1_seq, c2_seq, c3_seq)``
    :rtype: Tuple[bool, Sequence[float], Sequence[float], Optional[Sequence[float]]]
    :raises WrongDimensionErr: If sequences are empty or have mismatched lengths
    :raises InvalidCoordinateErr: If any coordinate is non-numeric or non-finite (NaN/Inf)
    """
    is_seq1 = isinstance(c1, (list, tuple))
    is_seq2 = isinstance(c2, (list, tuple))
    is_seq3 = isinstance(c3, (list, tuple)) if c3 is not None else is_seq1

    if is_seq1 != is_seq2 or (c3 is not None and is_seq1 != is_seq3):
        raise exceptions.WrongDimensionErr(
            "Coordinate arguments must all be scalars or all sequences of equal length."
        )

    is_scalar = not is_seq1
    if is_scalar:
        try:
            f1 = float(c1)  # type: ignore[arg-type]
            f2 = float(c2)  # type: ignore[arg-type]
            f3 = float(c3) if c3 is not None else None  # type: ignore[arg-type]
        except (ValueError, TypeError) as err:
            raise exceptions.InvalidCoordinateErr("Coordinate values must be numeric.") from err

        if not (math.isfinite(f1) and math.isfinite(f2) and (f3 is None or math.isfinite(f3))):
            raise exceptions.InvalidCoordinateErr(
                "Coordinate values must be finite numbers (NaN or Infinity is not allowed)."
            )
        return True, (f1,), (f2,), ((f3,) if f3 is not None else None)

    seq1: Sequence[float] = c1  # type: ignore[assignment]
    seq2: Sequence[float] = c2  # type: ignore[assignment]
    seq3: Optional[Sequence[float]] = c3  # type: ignore[assignment]

    n_len = len(seq1)
    if n_len == 0 or len(seq2) == 0 or (seq3 is not None and len(seq3) == 0):
        raise exceptions.WrongDimensionErr("Coordinate sequences cannot be empty.")

    if len(seq2) != n_len or (seq3 is not None and len(seq3) != n_len):
        raise exceptions.WrongDimensionErr("Coordinate sequences must have identical lengths.")

    try:
        norm1 = [float(x) for x in seq1]
        norm2 = [float(x) for x in seq2]
        norm3 = [float(x) for x in seq3] if seq3 is not None else None
    except (ValueError, TypeError) as err:
        raise exceptions.InvalidCoordinateErr(
            "Coordinate sequence elements must be numeric."
        ) from err

    if not (
        all(math.isfinite(x) for x in norm1)
        and all(math.isfinite(x) for x in norm2)
        and (norm3 is None or all(math.isfinite(x) for x in norm3))
    ):
        raise exceptions.InvalidCoordinateErr(
            "Coordinate sequence elements must be finite numbers (NaN or Infinity is not allowed)."
        )

    return False, norm1, norm2, norm3


class TransRO:
    """Coordinates geodetic transformations between Stereo 70 and ETRS89 coordinate systems."""

    def __init__(
        self,
        grid_filename: Optional[str] = None,
        telemetry_logger: Optional[Any] = None,
    ) -> None:
        """Initialize transformation engines with SPG binary grid and optional telemetry logger.

        :param grid_filename: Optional filename or path of the .spg grid file.
            If None, the grid is auto-discovered from pytransdatro/grids/.
        :type grid_filename: Optional[str]
        :param telemetry_logger: Optional asynchronous telemetry collector instance
        :type telemetry_logger: Optional[Any]
        """
        self._t_gr2d: Grid2D = Grid2D(grid_filename)
        self._t_h2d: Helmert2D = Helmert2D()
        self._p_st70: StereoProj = StereoProj()
        self._t_gr1d: Grid1D = Grid1D(grid_filename)
        self._collector: Optional[Any] = telemetry_logger

        self.grid2d_bounds_st70: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)
        self.grid2d_bounds_etrs89: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)
        self.grid1d_bounds_etrs89: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)
        self.grid1d_bounds_st70: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)
        self._init_dual_crs_extents()

    def with_collector(self, telemetry_logger: Optional[Any]) -> "TransRO":
        """Returns a shallow clone of TransRO with an alternative telemetry collector.

        Shares the pre-loaded 2D and 1D binary grids without memory duplication.
        """
        cloned = copy.copy(self)
        cloned._collector = telemetry_logger
        return cloned

    def _init_dual_crs_extents(self) -> None:
        """Computes and caches dual-CRS 6-point extended bounding boxes for Grid2D and Grid1D.

        Evaluates the 4 grid corners plus the 2 intersections with the Oblique Stereographic
        projection center (:math:`\\varphi_0 = 46^\\circ\\text{ N}, \\lambda_0 = 25^\\circ\\text{ E}`
        / :math:`E_0 = 500\\,000\\text{ m}`) to account for outward parallel curvature and
        meridian convergence.
        """
        g2 = self._t_gr2d
        g1 = self._t_gr1d
        e_center = self._p_st70.FALSE_E
        lon_center_rad = self._p_st70.ORIG_LON

        # 1. Grid2D: native Stereo70 meters -> 6-point extended ETRS89 radians
        self.grid2d_bounds_st70 = (g2.n_min, g2.n_max, g2.e_min, g2.e_max)
        g2_n6 = [g2.n_min, g2.n_min, g2.n_min, g2.n_max, g2.n_max, g2.n_max]
        g2_e6 = [g2.e_min, e_center, g2.e_max, g2.e_min, e_center, g2.e_max]
        g2_nh, g2_eh = self._t_h2d.trans(g2_n6, g2_e6, 1)
        g2_lats, g2_lons = self._p_st70.to_geo(g2_nh, g2_eh)
        self.grid2d_bounds_etrs89 = (
            min(g2_lats),
            max(g2_lats),
            min(g2_lons),
            max(g2_lons),
        )

        # 2. Grid1D: native ETRS89 radians -> 6-point extended Stereo70 meters
        lat1_min = math.radians(g1.n_min)
        lat1_max = math.radians(g1.n_max)
        lon1_min = math.radians(g1.e_min)
        lon1_max = math.radians(g1.e_max)
        self.grid1d_bounds_etrs89 = (lat1_min, lat1_max, lon1_min, lon1_max)

        g1_lat6 = [lat1_min, lat1_min, lat1_min, lat1_max, lat1_max, lat1_max]
        g1_lon6 = [lon1_min, lon_center_rad, lon1_max, lon1_min, lon_center_rad, lon1_max]
        g1_np, g1_ep = self._p_st70.to_grid(g1_lat6, g1_lon6)
        g1_ns, g1_es = self._t_h2d.trans(g1_np, g1_ep, -1)
        shift_margin_m = 10.0  # Exceeds max 2D grid shift (~4.4 m)
        self.grid1d_bounds_st70 = (
            min(g1_ns) - shift_margin_m,
            max(g1_ns) + shift_margin_m,
            min(g1_es) - shift_margin_m,
            max(g1_es) + shift_margin_m,
        )

    def _precheck_st70_point(self, n_val: float, e_val: float, is_3d: bool) -> None:
        """Pre-checks a single Stereo70 coordinate against Grid2D and (if 3D) Grid1D extents."""
        if not self._t_gr2d.covered_by_grid(n_val, e_val):
            raise exceptions.OutOfGrid2DErr(n_val, e_val, self._t_gr2d)
        if is_3d:
            n1_min, n1_max, e1_min, e1_max = self.grid1d_bounds_st70
            if not (n1_min <= n_val <= n1_max and e1_min <= e_val <= e1_max):
                raise exceptions.OutOfGrid1DErr(n_val, e_val, self._t_gr1d)

    def _precheck_etrs89_point(self, lat_rad: float, lon_rad: float, is_3d: bool) -> None:
        """Pre-checks a single ETRS89 radian coordinate against Grid2D and (if 3D) Grid1D extents."""
        lat2_min, lat2_max, lon2_min, lon2_max = self.grid2d_bounds_etrs89
        if not (lat2_min <= lat_rad <= lat2_max and lon2_min <= lon_rad <= lon2_max):
            raise exceptions.OutOfGrid2DErr(
                math.degrees(lat_rad), math.degrees(lon_rad), self._t_gr2d
            )
        if is_3d:
            lat1_min, lat1_max, lon1_min, lon1_max = self.grid1d_bounds_etrs89
            if not (lat1_min <= lat_rad <= lat1_max and lon1_min <= lon_rad <= lon1_max):
                raise exceptions.OutOfGrid1DErr(
                    math.degrees(lat_rad), math.degrees(lon_rad), self._t_gr1d
                )

    @overload
    def st70_to_etrs89(
        self,
        n: float,
        e: float,
        z: None = None,
        source: int = 0,
    ) -> Coord2D: ...

    @overload
    def st70_to_etrs89(
        self,
        n: float,
        e: float,
        z: float,
        source: int = 0,
    ) -> Coord3D: ...

    @overload
    def st70_to_etrs89(
        self,
        n: Sequence[float],
        e: Sequence[float],
        z: None = None,
        source: int = 0,
    ) -> Sequence2D: ...

    @overload
    def st70_to_etrs89(
        self,
        n: Sequence[float],
        e: Sequence[float],
        z: Sequence[float],
        source: int = 0,
    ) -> Sequence3D: ...

    def st70_to_etrs89(
        self,
        n: ScalarOrSeq,
        e: ScalarOrSeq,
        z: Optional[ScalarOrSeq] = None,
        source: int = 0,
    ) -> Union[Coord2D, Coord3D, Sequence2D, Sequence3D]:
        """Transforms Stereo70 grid coordinates to ETRS89 geographic coordinates.

        Supports both single coordinates (floats) and coordinate sequences (lists/tuples).
        :math:`(N, E[, H]) \\to (\\text{Lat}, \\text{Long}[, h])`

        :param n: Northing coordinate in meters (or sequence of northings)
        :type n: float or Sequence[float]
        :param e: Easting coordinate in meters (or sequence of eastings)
        :type e: float or Sequence[float]
        :param z: Normal elevation in meters (or sequence of elevations), optional
        :type z: float or Sequence[float], optional
        :param source: Integer source identifier (e.g. 0=unknown, 1=web_map, 2=rest_api)
        :type source: int, optional
        :return: (lat, lon[, h]) in ETRS89 coordinate reference system (radians, meters)
        :rtype: Tuple[float, float] or Tuple[float, float, float] or Sequence2D or Sequence3D
        :raises InvalidCoordinateErr: If any coordinate is non-numeric or non-finite (NaN/Inf)
        :raises WrongDimensionErr: If sequence inputs are empty or have mismatched lengths
        :raises OutOfGridErr: If (n, e) falls outside the active grid boundaries
        :raises NoDataGridErr: If subgrid interpolation encounters No Data (NaN or 999) values
        """
        is_scalar, n_seq, e_seq, z_seq = _validate_and_normalize_inputs(n, e, z)
        is_3d_bool = z_seq is not None

        for n_val, e_val in zip(n_seq, e_seq, strict=True):
            self._precheck_st70_point(n_val, e_val, is_3d_bool)

        r_n, r_e = self._t_gr2d.trans(n_seq, e_seq, -1)
        r_n, r_e = self._t_h2d.trans(r_n, r_e, 1)
        r_n, r_e = self._p_st70.to_geo(r_n, r_e)

        # If z value provided, apply 1D grid elevation correction: h = H + N
        r_z: Optional[list[float]] = None
        if z_seq is not None:
            _degrees = math.degrees
            n_deg = [_degrees(x) for x in r_n]
            e_deg = [_degrees(y) for y in r_e]
            (r_z,) = self._t_gr1d.trans(n_deg, e_deg, z_seq, 1)

        if self._collector is not None:
            try:
                is_3d = 1 if z_seq is not None else 0
                self._collector.record(1, r_n, r_e, is_3d=is_3d, source=source)
            except Exception as exc:
                logger.debug("Telemetry logging skipped during st70_to_etrs89: %s", exc)

        if r_z is not None:
            if is_scalar:
                return (r_n[0], r_e[0], r_z[0])
            return (r_n, r_e, r_z)
        else:
            if is_scalar:
                return (r_n[0], r_e[0])
            return (r_n, r_e)

    @overload
    def etrs89_to_st70(
        self,
        lat: float,
        lon: float,
        h: None = None,
        source: int = 0,
    ) -> Coord2D: ...

    @overload
    def etrs89_to_st70(
        self,
        lat: float,
        lon: float,
        h: float,
        source: int = 0,
    ) -> Coord3D: ...

    @overload
    def etrs89_to_st70(
        self,
        lat: Sequence[float],
        lon: Sequence[float],
        h: None = None,
        source: int = 0,
    ) -> Sequence2D: ...

    @overload
    def etrs89_to_st70(
        self,
        lat: Sequence[float],
        lon: Sequence[float],
        h: Sequence[float],
        source: int = 0,
    ) -> Sequence3D: ...

    def etrs89_to_st70(
        self,
        lat: ScalarOrSeq,
        lon: ScalarOrSeq,
        h: Optional[ScalarOrSeq] = None,
        source: int = 0,
    ) -> Union[Coord2D, Coord3D, Sequence2D, Sequence3D]:
        """Transforms ETRS89 geographic coordinates to Stereo70 grid coordinates.

        Supports both single coordinates (floats) and coordinate sequences (lists/tuples).
        :math:`(\\text{Lat}, \\text{Long}[, h]) \\to (N, E[, H])`

        :param lat: Latitude in radians (or sequence of latitudes)
        :type lat: float or Sequence[float]
        :param lon: Longitude in radians (or sequence of longitudes)
        :type lon: float or Sequence[float]
        :param h: Ellipsoidal elevation in meters (or sequence of elevations), optional
        :type h: float or Sequence[float], optional
        :param source: Integer source identifier (e.g. 0=unknown, 1=web_map, 2=rest_api)
        :type source: int, optional
        :return: (n, e[, z]) in Stereo70 coordinate reference system (meters)
        :rtype: Tuple[float, float] or Tuple[float, float, float] or Sequence2D or Sequence3D
        :raises InvalidCoordinateErr: If any coordinate is non-numeric or non-finite (NaN/Inf)
        :raises WrongDimensionErr: If sequence inputs are empty or have mismatched lengths
        :raises OutOfGridErr: If (lat, lon) falls outside the active grid boundaries
        :raises NoDataGridErr: If subgrid interpolation encounters No Data (NaN or 999) values
        """
        is_scalar, lat_seq, lon_seq, h_seq = _validate_and_normalize_inputs(lat, lon, h)
        is_3d_bool = h_seq is not None

        for lat_val, lon_val in zip(lat_seq, lon_seq, strict=True):
            self._precheck_etrs89_point(lat_val, lon_val, is_3d_bool)

        r_n, r_e = self._p_st70.to_grid(lat_seq, lon_seq)
        r_n, r_e = self._t_h2d.trans(r_n, r_e, -1)
        r_n, r_e = self._t_gr2d.trans(r_n, r_e, 1)

        # If h value provided, apply 1D grid elevation correction: H = h - N
        r_z: Optional[list[float]] = None
        if h_seq is not None:
            _degrees = math.degrees
            lat_deg = [_degrees(x) for x in lat_seq]
            lon_deg = [_degrees(y) for y in lon_seq]
            (r_z,) = self._t_gr1d.trans(lat_deg, lon_deg, h_seq, -1)

        if self._collector is not None:
            try:
                is_3d = 1 if h_seq is not None else 0
                self._collector.record(2, lat_seq, lon_seq, is_3d=is_3d, source=source)
            except Exception as exc:
                logger.debug("Telemetry logging skipped during etrs89_to_st70: %s", exc)

        if r_z is not None:
            if is_scalar:
                return (r_n[0], r_e[0], r_z[0])
            return (r_n, r_e, r_z)
        else:
            if is_scalar:
                return (r_n[0], r_e[0])
            return (r_n, r_e)

    def transform_batch(
        self,
        points: Sequence[Any],
        from_st70: bool = True,
        is_degrees: bool = False,
        source: int = 0,
    ) -> list[Tuple[list[float], Optional[exceptions.PyTransDatError]]]:
        """Executes fault-isolated batch coordinate transformation across a sequence of points.

        Implements the 4-stage pipeline:
            1. Per-point syntax, dimensionality (2D/3D), and 6-point Extended Bounding Box pre-check.
            2. Batch Oblique Stereographic projection and 2D Helmert transformation on survivors.
            3. Per-element fault-isolated Grid2D and Grid1D interpolation (catching curved-wedge
               ``OutOfGrid2DErr`` / ``OutOfGrid1DErr`` and ``NoDataGridErr`` without aborting valid points).
            4. Single batch telemetry aggregation and index-order result merging.

        :param points: Sequence of coordinate rows (each a 2D or 3D sequence of numbers)
        :type points: Sequence[Any]
        :param from_st70: True for Stereo70 -> ETRS89, False for ETRS89 -> Stereo70
        :type from_st70: bool
        :param is_degrees: True if ETRS89 coordinates are in decimal degrees, False for radians
        :type is_degrees: bool
        :param source: Telemetry client source ID
        :type source: int
        :return: List of ``(output_or_echo_coords, domain_error_or_none)`` in original input order
        :rtype: list[Tuple[list[float], Optional[exceptions.PyTransDatError]]]
        """
        n_total = len(points)
        results: list[Optional[Tuple[list[float], Optional[exceptions.PyTransDatError]]]] = [
            None
        ] * n_total

        surv_idxs: list[int] = []
        surv_c1: list[float] = []
        surv_c2: list[float] = []
        surv_c3: list[Optional[float]] = []
        surv_orig: list[list[float]] = []

        # Stage 1: Syntax, Dimension, and 6-Point Extended Bounding Box Pre-Check
        for idx, raw_pt in enumerate(points):
            if not isinstance(raw_pt, (list, tuple)):
                results[idx] = ([], exceptions.InvalidCoordinateErr("Invalid coordinate data"))
                continue

            coord_list: list[float] = []
            invalid_elem = False
            for elem in raw_pt:
                if isinstance(elem, bool) or elem is None:
                    invalid_elem = True
                    break
                try:
                    val = float(elem)
                except (ValueError, TypeError):
                    invalid_elem = True
                    break
                if not math.isfinite(val):
                    invalid_elem = True
                    break
                coord_list.append(val)

            if invalid_elem:
                results[idx] = ([], exceptions.InvalidCoordinateErr("Invalid coordinate data"))
                continue

            dim = len(coord_list)
            if dim < 2 or dim > 3:
                results[idx] = (coord_list, exceptions.WrongDimensionErr("Wrong dimension"))
                continue

            c1, c2 = coord_list[0], coord_list[1]
            c3 = coord_list[2] if dim == 3 else None
            is_3d = c3 is not None

            try:
                if from_st70:
                    self._precheck_st70_point(c1, c2, is_3d)
                    surv_c1.append(c1)
                    surv_c2.append(c2)
                else:
                    lat_rad = math.radians(c1) if is_degrees else c1
                    lon_rad = math.radians(c2) if is_degrees else c2
                    self._precheck_etrs89_point(lat_rad, lon_rad, is_3d)
                    surv_c1.append(lat_rad)
                    surv_c2.append(lon_rad)
            except exceptions.OutOfGridErr as err:
                results[idx] = (coord_list, err)
                continue

            surv_idxs.append(idx)
            surv_c3.append(c3)
            surv_orig.append(coord_list)

        if not surv_idxs:
            return [r for r in results if r is not None]

        # Track valid points for single-call batch telemetry
        telem_2d_lats: list[float] = []
        telem_2d_lons: list[float] = []
        telem_3d_lats: list[float] = []
        telem_3d_lons: list[float] = []

        # Stages 2 & 3: Vectorized Projection/Helmert + Fault-Isolated Grid Interpolation
        if from_st70:
            s2_idxs: list[int] = []
            s2_rn: list[float] = []
            s2_re: list[float] = []
            s2_c3: list[Optional[float]] = []
            s2_orig: list[list[float]] = []

            for k, orig_idx in enumerate(surv_idxs):
                try:
                    rn_1, re_1 = self._t_gr2d.trans([surv_c1[k]], [surv_c2[k]], -1)
                    s2_idxs.append(orig_idx)
                    s2_rn.append(rn_1[0])
                    s2_re.append(re_1[0])
                    s2_c3.append(surv_c3[k])
                    s2_orig.append(surv_orig[k])
                except (exceptions.OutOfGridErr, exceptions.NoDataGridErr) as err:
                    results[orig_idx] = (surv_orig[k], err)

            if s2_idxs:
                h_n, h_e = self._t_h2d.trans(s2_rn, s2_re, 1)
                geo_lats, geo_lons = self._p_st70.to_geo(h_n, h_e)

                for m, orig_idx in enumerate(s2_idxs):
                    lat_rad = geo_lats[m]
                    lon_rad = geo_lons[m]
                    z_val = s2_c3[m]
                    lat_out = math.degrees(lat_rad) if is_degrees else lat_rad
                    lon_out = math.degrees(lon_rad) if is_degrees else lon_rad

                    if z_val is not None:
                        try:
                            (rz_1,) = self._t_gr1d.trans(
                                [math.degrees(lat_rad)],
                                [math.degrees(lon_rad)],
                                [z_val],
                                1,
                            )
                        except (exceptions.OutOfGridErr, exceptions.NoDataGridErr) as err:
                            results[orig_idx] = (s2_orig[m], err)
                            continue

                        results[orig_idx] = ([lat_out, lon_out, rz_1[0]], None)
                        telem_3d_lats.append(lat_rad)
                        telem_3d_lons.append(lon_rad)
                    else:
                        results[orig_idx] = ([lat_out, lon_out], None)
                        telem_2d_lats.append(lat_rad)
                        telem_2d_lons.append(lon_rad)

            op_code = 1
        else:
            proj_n, proj_e = self._p_st70.to_grid(surv_c1, surv_c2)
            helm_n, helm_e = self._t_h2d.trans(proj_n, proj_e, -1)

            for k, orig_idx in enumerate(surv_idxs):
                lat_rad = surv_c1[k]
                lon_rad = surv_c2[k]
                h_val = surv_c3[k]

                try:
                    rn_1, re_1 = self._t_gr2d.trans([helm_n[k]], [helm_e[k]], 1)
                except (exceptions.OutOfGridErr, exceptions.NoDataGridErr) as err:
                    results[orig_idx] = (surv_orig[k], err)
                    continue

                if h_val is not None:
                    try:
                        (rz_1,) = self._t_gr1d.trans(
                            [math.degrees(lat_rad)],
                            [math.degrees(lon_rad)],
                            [h_val],
                            -1,
                        )
                    except (exceptions.OutOfGridErr, exceptions.NoDataGridErr) as err:
                        results[orig_idx] = (surv_orig[k], err)
                        continue

                    results[orig_idx] = ([rn_1[0], re_1[0], rz_1[0]], None)
                    telem_3d_lats.append(lat_rad)
                    telem_3d_lons.append(lon_rad)
                else:
                    results[orig_idx] = ([rn_1[0], re_1[0]], None)
                    telem_2d_lats.append(lat_rad)
                    telem_2d_lons.append(lon_rad)

            op_code = 2

        # Stage 4: Single Batch Telemetry Call
        if self._collector is not None:
            try:
                if telem_2d_lats:
                    self._collector.record(
                        op_code, telem_2d_lats, telem_2d_lons, is_3d=0, source=source
                    )
                if telem_3d_lats:
                    self._collector.record(
                        op_code, telem_3d_lats, telem_3d_lons, is_3d=1, source=source
                    )
            except Exception as exc:
                logger.debug("Telemetry logging skipped during transform_batch: %s", exc)

        return [r for r in results if r is not None]
