#!/usr/bin/env python3
"""Draw the ViroSync workflow map that the README shows.

The map follows the data flow of the single-genome flow in
``virosync.orchestration._flows.single_genome``: a step sits after the steps
whose output it reads. Before it writes the PNG, the script checks the drawing.
Every station must touch exactly the lines it lists, every line must attach
where it says it does, no two stations or lines may touch elsewhere, and every
label must sit inside its phase band, nearer its own station than any other and
clear of all other marks. All text must render at 12 px or more in the README.
The map needs the Liberation Sans font (Debian and Ubuntu package
``fonts-liberation``).
"""

from __future__ import annotations

import argparse
import functools
import io
import itertools
import logging
import math
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, TypeAlias

import numpy as np
from matplotlib import font_manager, rc_context
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import FancyBboxPatch, PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.text import Text
from matplotlib.transforms import Bbox
from PIL import Image

logger = logging.getLogger(__name__)

DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "docs" / "virosync_workflow.png"
README_WIDTH_PX = 830  # README column of GitHub in a desktop browser
MIN_TEXT_PX = 12.0
WIDTH_PT = 672.0
HEIGHT_PT = 228.0
DPI = 216  # three pixels per point, so whole-point sizes give whole pixels
FONT = "Liberation Sans"
LABEL_PT = 9.75
HEADING_PT = 10.0
TEXT_COLOR = "#1F2328"
MUTED_COLOR = "#57606A"
BAND_COLOR = "#F2F4F7"
OUTLINE_COLOR = "#24292F"

STOP_RADIUS = 4.4
STOP_EDGE = 1.7
STOP_HALF = STOP_RADIUS + STOP_EDGE / 2
LABEL_GAP = 4.0
CORNER_RADIUS = 7.0
CLEARANCE = 1.5  # free space around every label, stop and line
ROW_GAP = 4.0  # horizontal gap between labels that share a row
EDGE_MARGIN = 4.0  # distance between a vertical line and a band edge
JUNCTION_PT = 15.0  # a branch may run along its parent line this far from its end
ON_LINE_PT = 0.5  # how far an end may sit from the stop or line it attaches to
SAMPLE_PT = 0.5  # spacing of the points that stand for a line in the checks

BAND_BOTTOM = 24.0
BAND_TOP = HEIGHT_PT - 4.0
Y_MAIN = 83.0
Y_HOST = Y_MAIN - 11.0
Y_OPTIONAL = Y_MAIN + 50.0
Y_ANNOTATION = Y_OPTIONAL + 18.0
Y_LEGEND = 12.0
ANNOTATION_RISER_X = 462.0
OPTIONAL_RISER_X = 482.0

Point: TypeAlias = tuple[float, float]


@dataclass(frozen=True, slots=True)
class Track:
    """One line of the map: its name, legend text, colour, width and dashes."""

    name: str
    legend: str
    color: str
    width: float
    dashes: tuple[float, float] | None = None


MAIN = Track("main", "Main path", "#24292F", 4.6)
HOST = Track("host", "Host and taxonomy evidence", "#009E73", 4.0)
OPTIONAL = Track("optional", "Optional, off by default", "#CC79A7", 3.2, (2.2, 1.3))
ANNOTATION = Track("annotation", "Annotation only, no score effect", "#D55E00", 2.2)


@dataclass(frozen=True, slots=True)
class Phase:
    """A phase band: its number, its title and its x range."""

    number: int
    title: str
    x0: float
    x1: float

    @property
    def box(self) -> Bbox:
        """Return the band rectangle."""
        return Bbox.from_extents(self.x0, BAND_BOTTOM, self.x1, BAND_TOP)


PHASES = (
    Phase(0, "Preprocessing", 8.0, 76.0),
    Phase(1, "Markers and seeds", 79.0, 258.0),
    Phase(2, "Boundaries", 261.0, 492.0),
    Phase(3, "Scoring and output", 495.0, 664.0),
)


@dataclass(frozen=True, slots=True)
class Stop:
    """A station: its label, its position and the lines it serves.

    A ``dot`` sits on one line, an ``interchange`` spans the lines between
    ``y_low`` and ``y_high`` and a ``terminal`` starts or ends the main path.
    Every kind is drawn as a rounded bar of radius STOP_HALF around the vertical
    segment from ``y_low`` to ``y_high``. An optional stop is off in a default run.
    """

    label: str
    x: float
    y_low: float
    y_high: float
    tracks: tuple[Track, ...]
    side: Literal["above", "below"]
    kind: Literal["dot", "interchange", "terminal"] = "dot"
    align: Literal["left", "center", "right"] = "center"
    is_optional: bool = False


GENOME = Stop("Genome\nassembly", 20.0, Y_MAIN, Y_OPTIONAL, (MAIN, OPTIONAL), "below", "terminal", "left")
VALIDATION = Stop(
    "Marker\nvalidation\nPfam, DIAMOND", 146.0, Y_HOST, Y_OPTIONAL, (MAIN, HOST, OPTIONAL), "below", "interchange"
)
SYNTHESIS = Stop(
    "Evidence\nsynthesis,\nconfidence tiers", 557.0, Y_HOST, Y_OPTIONAL, (MAIN, HOST, OPTIONAL), "below", "interchange"
)
OUTPUTS = Stop(
    "BED, GFF3,\nTSV, FASTA,\nreport", 653.0, Y_MAIN, Y_ANNOTATION, (MAIN, ANNOTATION), "below", "terminal", "right"
)

STOPS = (
    GENOME,
    Stop("Gene\nprediction", 51.0, Y_MAIN, Y_MAIN, (MAIN,), "above"),
    Stop("Codon-aware\nDNA scan", 112.0, Y_OPTIONAL, Y_OPTIONAL, (OPTIONAL,), "above"),
    Stop("Marker\nHMM scan", 112.0, Y_MAIN, Y_MAIN, (MAIN,), "above"),
    VALIDATION,
    Stop("Host\nsignature\nmodel", 210.0, Y_HOST, Y_HOST, (HOST,), "below"),
    Stop("Candidate\nregions", 229.0, Y_MAIN, Y_MAIN, (MAIN,), "above"),
    Stop("Host\nsignature\ntrim", 288.0, Y_HOST, Y_MAIN, (MAIN, HOST), "above", "interchange"),
    Stop("Extension,\ngene taxonomy,\nDIAMOND", 324.0, Y_HOST, Y_MAIN, (MAIN, HOST), "below", "interchange"),
    Stop("Boundary\nrefinement", 358.0, Y_HOST, Y_MAIN, (MAIN, HOST), "above", "interchange"),
    Stop("Gap merge,\nTIR ends, splits", 396.0, Y_HOST, Y_MAIN, (MAIN, HOST), "below", "interchange"),
    Stop("Composition\nGC, k-mer", 430.0, Y_MAIN, Y_MAIN, (MAIN,), "above"),
    Stop("Direct\nrepeats", 472.0, Y_ANNOTATION, Y_ANNOTATION, (ANNOTATION,), "above"),
    Stop("Domains,\nstructures,\nphylogeny", 506.0, Y_OPTIONAL, Y_OPTIONAL, (OPTIONAL,), "below", "dot", "left"),
    SYNTHESIS,
    Stop("Integration\ngenes", 570.0, Y_ANNOTATION, Y_ANNOTATION, (ANNOTATION,), "above", "dot", "right"),
    Stop("Acceptance,\nANI class\ntransfer", 594.0, Y_MAIN, Y_MAIN, (MAIN,), "above"),
    Stop(
        "Gene-function\nlabels, bundle\n≥ v1.2.0",
        622.0,
        Y_ANNOTATION,
        Y_ANNOTATION,
        (ANNOTATION,),
        "above",
        is_optional=True,
    ),
    OUTPUTS,
)


@dataclass(frozen=True, slots=True)
class Route:
    """One drawn line: its track, its points and what each end attaches to.

    An end attaches to a stop when it lies on that stop, or to a track when the
    line branches from that track's line there. None leaves an end free, as for
    the legend handles.
    """

    track: Track
    points: tuple[Point, ...]
    start: Stop | Track | None
    end: Stop | Track | None


def branch(track: Track, riser_x: float, y: float, end: Stop) -> Route:
    """Return a line that leaves the main path at *riser_x*, rises to *y* and runs to *end*."""
    points = ((riser_x - CORNER_RADIUS, Y_MAIN), (riser_x, Y_MAIN), (riser_x, y), (end.x, y))
    return Route(track, points, MAIN, end)


# Each line starts on the stop or line that feeds it and ends where it is read.
ROUTES = (
    Route(OPTIONAL, ((GENOME.x + STOP_HALF, Y_OPTIONAL), (VALIDATION.x, Y_OPTIONAL)), GENOME, VALIDATION),
    branch(OPTIONAL, OPTIONAL_RISER_X, Y_OPTIONAL, SYNTHESIS),
    Route(HOST, ((VALIDATION.x, Y_HOST), (SYNTHESIS.x, Y_HOST)), VALIDATION, SYNTHESIS),
    branch(ANNOTATION, ANNOTATION_RISER_X, Y_ANNOTATION, OUTPUTS),
    Route(MAIN, ((GENOME.x, Y_MAIN), (OUTPUTS.x, Y_MAIN)), GENOME, OUTPUTS),
)


@dataclass(frozen=True, slots=True)
class Label:
    """A drawn text, the stop it names and the band it must stay in."""

    text: Text
    stop: Stop | None
    phase: Phase | None

    @property
    def name(self) -> str:
        """Return the text on one line."""
        return self.text.get_text().replace("\n", " ")


@dataclass(frozen=True, slots=True)
class Drawn:
    """Everything the layout check reads: labels, stops and routes."""

    labels: list[Label]
    stops: list[Stop]
    routes: list[Route]


def draw_map(stops: Sequence[Stop] = STOPS, routes: Sequence[Route] = ROUTES) -> tuple[Figure, Drawn]:
    """Draw the bands, routes, stops, labels and legend on a new figure."""
    figure = Figure(figsize=(WIDTH_PT / 72, HEIGHT_PT / 72), dpi=DPI, facecolor="white")
    FigureCanvasAgg(figure)
    axes = figure.add_axes((0, 0, 1, 1))
    axes.set_xlim(0, WIDTH_PT)
    axes.set_ylim(0, HEIGHT_PT)
    axes.set_axis_off()
    labels = draw_phases(axes)
    for route in routes:
        draw_route(axes, route)
    for stop in stops:
        draw_stop(axes, stop)
        labels.append(Label(draw_label(axes, stop, routes), stop, phase_of(stop.x)))
    legend_labels, legend_stop, legend_routes = draw_legend(axes)
    labels += legend_labels
    return figure, Drawn(labels, [*stops, legend_stop], [*routes, *legend_routes])


def draw_phases(axes: Axes) -> list[Label]:
    """Draw one light band per phase with its number and title on top."""
    headings = []
    for phase in PHASES:
        band = phase.box
        axes.add_patch(
            FancyBboxPatch(
                (band.x0, band.y0),
                band.width,
                band.height,
                boxstyle="round,pad=0,rounding_size=4",
                facecolor=BAND_COLOR,
                edgecolor="none",
                zorder=0,
            )
        )
        centre = (phase.x0 + phase.x1) / 2
        number = write(axes, centre, BAND_TOP - 5.0, f"Phase {phase.number}", weight="bold", size=HEADING_PT)
        title = write(axes, centre, BAND_TOP - 17.5, phase.title, color=MUTED_COLOR)
        headings += [Label(number, None, phase), Label(title, None, phase)]
    return headings


def draw_route(axes: Axes, route: Route) -> None:
    """Draw one line with rounded bends."""
    patch = PathPatch(
        rounded_path(route.points),
        facecolor="none",
        edgecolor=route.track.color,
        linewidth=route.track.width,
        capstyle="butt",
        joinstyle="round",
        zorder=3 if route.track == MAIN else 2,
    )
    if route.track.dashes is not None:
        patch.set_linestyle((0, route.track.dashes))
    axes.add_patch(patch)


def rounded_path(points: Sequence[Point]) -> MplPath:
    """Join *points* with straight segments and round each bend.

    A bend takes up to CORNER_RADIUS of each segment it touches, all of an end
    segment and half of a segment it shares with another bend.
    """
    vertices = [points[0]]
    codes = [MplPath.MOVETO]
    last = len(points) - 1
    for index in range(1, last):
        before, corner, after = points[index - 1], points[index], points[index + 1]
        reach = min(
            CORNER_RADIUS,
            math.dist(before, corner) / (1 if index == 1 else 2),
            math.dist(corner, after) / (1 if index == last - 1 else 2),
        )
        vertices += [toward(corner, before, reach), corner, toward(corner, after, reach)]
        codes += [MplPath.LINETO, MplPath.CURVE3, MplPath.CURVE3]
    vertices.append(points[-1])
    codes.append(MplPath.LINETO)
    return MplPath(vertices, codes)


def toward(start: Point, end: Point, distance: float) -> Point:
    """Return the point *distance* from *start* toward *end*."""
    fraction = distance / math.dist(start, end)
    return (start[0] + (end[0] - start[0]) * fraction, start[1] + (end[1] - start[1]) * fraction)


def draw_stop(axes: Axes, stop: Stop) -> None:
    """Draw one stop: a ring on one line, a capsule across lines or a terminal bar."""
    fill = OUTLINE_COLOR if stop.kind == "terminal" else "white"
    edge = OUTLINE_COLOR
    if stop.kind == "dot":
        edge = OPTIONAL.color if stop.is_optional else stop.tracks[0].color
    axes.add_patch(
        FancyBboxPatch(
            (stop.x - STOP_RADIUS, stop.y_low - STOP_RADIUS),
            2 * STOP_RADIUS,
            stop.y_high - stop.y_low + 2 * STOP_RADIUS,
            boxstyle=f"round,pad=0,rounding_size={STOP_RADIUS}",
            facecolor=fill,
            edgecolor=edge,
            linewidth=STOP_EDGE,
            zorder=4,
        )
    )


def stop_box(stop: Stop) -> Bbox:
    """Return the rectangle around a stop, outline included."""
    return Bbox.from_extents(stop.x - STOP_HALF, stop.y_low - STOP_HALF, stop.x + STOP_HALF, stop.y_high + STOP_HALF)


def draw_label(axes: Axes, stop: Stop, routes: Sequence[Route]) -> Text:
    """Write the label of *stop* above it, or below the lowest line under it."""
    offset = STOP_HALF + LABEL_GAP
    x = {"left": stop.x - STOP_HALF, "center": stop.x, "right": stop.x + STOP_HALF}[stop.align]
    if stop.side == "above":
        return write(axes, x, stop.y_high + offset, stop.label, ha=stop.align, va="bottom")
    floor = min([stop.y_low, *lines_below(stop, routes)]) if stop.y_low <= Y_MAIN else stop.y_low
    return write(axes, x, floor - offset, stop.label, ha=stop.align, va="top")


def lines_below(stop: Stop, routes: Sequence[Route]) -> Iterator[float]:
    """Yield the y of every horizontal route segment that passes under *stop*."""
    for route in routes:
        for (x0, y0), (x1, y1) in itertools.pairwise(route.points):
            if y0 == y1 < stop.y_low and min(x0, x1) <= stop.x <= max(x0, x1):
                yield y0


def draw_legend(axes: Axes) -> tuple[list[Label], Stop, list[Route]]:
    """Draw the key in one row under the bands: four lines and the shared step."""
    handle = 22.0
    x = 14.0
    labels, routes = [], []
    for track in (MAIN, HOST, OPTIONAL, ANNOTATION):
        routes.append(Route(track, ((x, Y_LEGEND), (x + handle, Y_LEGEND)), None, None))
        text = write(axes, x + handle + 5.0, Y_LEGEND, track.legend, ha="left", va="center_baseline")
        labels.append(Label(text, None, None))
        x = text_box(axes.figure, text).x1 + 14.0
    sample = Stop("", x + 13.0, Y_LEGEND - 3.6, Y_LEGEND + 3.6, (MAIN, HOST), "above", "interchange")
    routes += [
        Route(MAIN, ((x, sample.y_high), (x + 26.0, sample.y_high)), None, None),
        Route(HOST, ((x, sample.y_low), (x + 26.0, sample.y_low)), None, None),
    ]
    for route in routes:
        draw_route(axes, route)
    draw_stop(axes, sample)
    text = write(axes, x + 31.0, Y_LEGEND, "Shared step", ha="left", va="center_baseline")
    labels.append(Label(text, None, None))
    return labels, sample, routes


def text_box(figure: Figure, text: Text) -> Bbox:
    """Return the box *text* covers, in data units."""
    renderer = figure.canvas.get_renderer()
    return text.get_window_extent(renderer).transformed(figure.axes[0].transData.inverted())


def write(
    axes: Axes,
    x: float,
    y: float,
    content: str,
    *,
    ha: str = "center",
    va: str = "top",
    color: str = TEXT_COLOR,
    size: float = LABEL_PT,
    weight: str = "normal",
) -> Text:
    """Place one text in the house font."""
    return axes.text(
        x,
        y,
        content,
        ha=ha,
        va=va,
        color=color,
        fontsize=size,
        fontweight=weight,
        fontfamily=font_family(),
        linespacing=1.15,
        multialignment=ha,
        zorder=5,
    )


@functools.cache
def font_family() -> str:
    """Return FONT, or raise when it is missing, so the layout never shifts silently."""
    for weight in ("normal", "bold"):
        properties = font_manager.FontProperties(family=FONT, weight=weight)
        try:
            font_manager.findfont(properties, fallback_to_default=False)
        except ValueError as err:
            raise FileNotFoundError(f"the map needs the {FONT} font (Debian and Ubuntu: fonts-liberation)") from err
    return FONT


def phase_of(x: float) -> Phase | None:
    """Return the phase band that contains *x*, or None between bands."""
    return next((phase for phase in PHASES if phase.x0 <= x <= phase.x1), None)


def layout_problems(figure: Figure, drawn: Drawn) -> list[str]:
    """List every place where the drawing breaks the map's rules."""
    return topology_problems(drawn.stops, drawn.routes) + label_problems(figure, drawn)


def topology_problems(stops: Sequence[Stop], routes: Sequence[Route]) -> list[str]:
    """List stops and lines that touch where the map says they do not, or miss where it says they do."""
    sampled = [(route, line_points(route)) for route in routes]
    problems = []
    for stop in stops:
        problems += connection_problems(stop, sampled)
        problems += placement_problems(stop_name(stop), stop_box(stop))
    for first, second in itertools.combinations(stops, 2):
        if stops_gap(first, second) < CLEARANCE:
            problems.append(f"stops {stop_name(first)!r} and {stop_name(second)!r} touch")
    for route, _points in sampled:
        problems += attachment_problems(route, sampled)
    for (first, points_a), (second, points_b) in itertools.combinations(sampled, 2):
        if lines_touch(first, points_a, second, points_b):
            problems.append(f"the {first.track.name} line and the {second.track.name} line touch")
    return problems + riser_problems(routes)


def connection_problems(stop: Stop, sampled: Sequence[tuple[Route, np.ndarray]]) -> list[str]:
    """Compare the lines that pass through *stop* with the lines it lists."""
    problems = [
        f"{stop_name(stop)!r} misses the {track.name} line"
        for track in stop.tracks
        if not any(route.track == track and point_gap(stop, points) <= ON_LINE_PT for route, points in sampled)
    ]
    for route, points in sampled:
        if route.track in stop.tracks:
            continue
        gap = point_gap(stop, points)
        if gap <= ON_LINE_PT:
            problems.append(f"the {route.track.name} line runs through {stop_name(stop)!r}")
        elif gap < route.track.width / 2 + CLEARANCE:
            problems.append(f"the {route.track.name} line grazes {stop_name(stop)!r}")
    return problems


def attachment_problems(route: Route, sampled: Sequence[tuple[Route, np.ndarray]]) -> list[str]:
    """Return a problem for each end of *route* that is not where the route says it attaches."""
    problems = []
    for anchor, (x, y) in ((route.start, route.points[0]), (route.end, route.points[-1])):
        if isinstance(anchor, Stop) and point_gap(anchor, np.array([[x, y]])) > ON_LINE_PT:
            problems.append(f"the {route.track.name} line is not attached to {stop_name(anchor)!r}")
        if isinstance(anchor, Track) and not any(
            other.track == anchor and other is not route and nearest(points, x, y) <= ON_LINE_PT
            for other, points in sampled
        ):
            problems.append(f"the {route.track.name} line is not attached to the {anchor.name} line")
    return problems


def placement_problems(name: str, box: Bbox) -> list[str]:
    """Return a problem when *box* straddles a band edge or leaves the canvas."""
    canvas = Bbox.from_extents(2.0, 2.0, WIDTH_PT - 2.0, HEIGHT_PT - 2.0)
    if not inside(box, canvas, 0.0):
        return [f"{name!r} leaves the canvas"]
    if any(inside(box, phase.box, CLEARANCE) for phase in PHASES):
        return []
    if any(padded(box).overlaps(phase.box) for phase in PHASES):
        return [f"{name!r} straddles a band edge"]
    return []


def riser_problems(routes: Sequence[Route]) -> list[str]:
    """Return a problem for each vertical line inside the bands that runs along a band edge."""
    return [
        f"the {route.track.name} line runs along the edge of phase {phase.number}"
        for route in routes
        for (x0, y0), (x1, y1) in itertools.pairwise(route.points)
        if x0 == x1 and max(y0, y1) > BAND_BOTTOM and min(y0, y1) < BAND_TOP
        for phase in PHASES
        if min(abs(x0 - phase.x0), abs(x0 - phase.x1)) < EDGE_MARGIN + route.track.width / 2
    ]


def label_problems(figure: Figure, drawn: Drawn) -> list[str]:
    """List labels that are too small, misplaced, or touch another mark."""
    sampled = [(route, line_points(route)) for route in drawn.routes]
    boxes = [(label, text_box(figure, label.text)) for label in drawn.labels]
    problems = [problem for label, _ in boxes if (problem := small_text(label))]
    for label, box in boxes:
        problems += placement_problems(label.name, box)
        if label.phase is not None and not inside(box, label.phase.box, CLEARANCE):
            problems.append(f"{label.name!r} leaves phase {label.phase.number}")
        problems += [
            f"{label.name!r} touches the {route.track.name} line"
            for route, points in sampled
            if box_reaches(box, points, route.track.width / 2 + CLEARANCE)
        ]
        problems += [
            f"{label.name!r} touches {stop_name(stop)!r}" for stop in drawn.stops if box_gap(stop, box) < CLEARANCE
        ]
        if label.stop is not None:
            problems += nearest_stop_problems(label, box, drawn.stops)
    for (first, box_a), (second, box_b) in itertools.combinations(boxes, 2):
        if too_close(box_a, box_b):
            problems.append(f"{first.name!r} and {second.name!r} are too close")
    return problems


def small_text(label: Label) -> str:
    """Return a problem when the label renders below the minimum size, else ''."""
    pixels = label.text.get_fontsize() * README_WIDTH_PX / WIDTH_PT
    if pixels >= MIN_TEXT_PX:
        return ""
    return f"{label.name!r} renders at {pixels:.1f} px"


def nearest_stop_problems(label: Label, box: Bbox, stops: Sequence[Stop]) -> list[str]:
    """Return a problem when another stop is nearer to the label than its own."""
    own = box_gap(label.stop, box)
    return [
        f"{label.name!r} is nearer to {stop_name(stop)!r} than to its own stop"
        for stop in stops
        if stop is not label.stop and box_gap(stop, box) < own
    ]


def too_close(first: Bbox, second: Bbox) -> bool:
    """Return whether two label boxes touch, or share a row without a clear gap."""
    if padded(first).overlaps(second):
        return True
    same_row = first.y0 < second.y1 and second.y0 < first.y1
    gap = max(second.x0 - first.x1, first.x0 - second.x1)
    return same_row and gap < ROW_GAP


def stop_name(stop: Stop) -> str:
    """Return the label of *stop* on one line."""
    return stop.label.replace("\n", " ") or "legend stop"


def point_gap(stop: Stop, points: np.ndarray) -> float:
    """Return the smallest distance from *points* to the outline of *stop*; negative inside."""
    dy = np.maximum(np.maximum(stop.y_low - points[:, 1], points[:, 1] - stop.y_high), 0.0)
    return float((np.hypot(points[:, 0] - stop.x, dy) - STOP_HALF).min())


def box_gap(stop: Stop, box: Bbox) -> float:
    """Return the distance from *box* to the outline of *stop*; negative when they overlap."""
    dx = max(box.x0 - stop.x, stop.x - box.x1, 0.0)
    dy = max(box.y0 - stop.y_high, stop.y_low - box.y1, 0.0)
    return math.hypot(dx, dy) - STOP_HALF


def stops_gap(first: Stop, second: Stop) -> float:
    """Return the distance between the outlines of two stops; negative when they overlap."""
    dy = max(second.y_low - first.y_high, first.y_low - second.y_high, 0.0)
    return math.hypot(first.x - second.x, dy) - 2 * STOP_HALF


def nearest(points: np.ndarray, x: float, y: float) -> float:
    """Return the distance from (x, y) to the nearest of *points*."""
    return float(np.hypot(points[:, 0] - x, points[:, 1] - y).min())


def inside(box: Bbox, frame: Bbox, margin: float) -> bool:
    """Return whether *box* lies inside *frame* with *margin* to spare."""
    return (
        box.x0 >= frame.x0 + margin
        and box.x1 <= frame.x1 - margin
        and box.y0 >= frame.y0 + margin
        and box.y1 <= frame.y1 - margin
    )


def padded(box: Bbox, margin: float = CLEARANCE) -> Bbox:
    """Grow *box* by *margin* on every side."""
    return Bbox.from_extents(box.x0 - margin, box.y0 - margin, box.x1 + margin, box.y1 + margin)


def box_reaches(box: Bbox, points: np.ndarray, reach: float) -> bool:
    """Return whether any of *points* lies within *reach* of *box*."""
    grown = padded(box, reach)
    xs, ys = points[:, 0], points[:, 1]
    return bool(((xs >= grown.x0) & (xs <= grown.x1) & (ys >= grown.y0) & (ys <= grown.y1)).any())


def lines_touch(first: Route, first_points: np.ndarray, second: Route, second_points: np.ndarray) -> bool:
    """Return whether two lines come within the clearance of each other outside their declared junctions."""
    points_a = outside_junctions(first, first_points, second.track)
    points_b = outside_junctions(second, second_points, first.track)
    if not len(points_a) or not len(points_b):
        return False
    gaps = np.hypot(points_a[:, None, 0] - points_b[None, :, 0], points_a[:, None, 1] - points_b[None, :, 1])
    return bool(gaps.min() < (first.track.width + second.track.width) / 2 + CLEARANCE)


def outside_junctions(route: Route, points: np.ndarray, other: Track) -> np.ndarray:
    """Drop the points of *route* near an end where it branches from the line of *other*."""
    keep = np.ones(len(points), dtype=bool)
    for anchor, (x, y) in ((route.start, route.points[0]), (route.end, route.points[-1])):
        if isinstance(anchor, Track) and anchor == other:
            keep &= np.hypot(points[:, 0] - x, points[:, 1] - y) > JUNCTION_PT
    return points[keep]


def line_points(route: Route) -> np.ndarray:
    """Return points every SAMPLE_PT along the centre line of *route*, curves included."""
    path = rounded_path(route.points)
    return np.array(list(sample_path(path, SAMPLE_PT)))


def sample_path(path: MplPath, step: float) -> Iterator[Point]:
    """Yield points every *step* along *path*, with curves cut into segments."""
    for polygon in path.to_polygons(closed_only=False):
        for (x0, y0), (x1, y1) in itertools.pairwise(polygon):
            count = max(1, math.ceil(math.hypot(x1 - x0, y1 - y0) / step))
            for index in range(count + 1):
                yield (x0 + (x1 - x0) * index / count, y0 + (y1 - y0) * index / count)


def planned_size() -> tuple[int, int]:
    """Return the PNG size in pixels that the canvas and DPI call for."""
    return (round(WIDTH_PT / 72 * DPI), round(HEIGHT_PT / 72 * DPI))


def render_png(figure: Figure) -> bytes:
    """Render the figure to PNG bytes, ignoring any savefig.bbox setting."""
    buffer = io.BytesIO()
    with rc_context({"savefig.bbox": "standard"}):
        figure.savefig(buffer, dpi=DPI, format="png", facecolor="white", metadata={"Software": None})
    return buffer.getvalue()


def size_problems(png: bytes) -> list[str]:
    """Return a problem when the PNG is not the planned size."""
    with Image.open(io.BytesIO(png)) as image:
        width, height = image.size
    expected = planned_size()
    if (width, height) == expected:
        return []
    return [f"the PNG is {width} x {height} px, not {expected[0]} x {expected[1]}"]


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    """Read the output path."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="PNG to write")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Draw the map, check it, and write the PNG only when every check passes."""
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    figure, drawn = draw_map()
    problems = layout_problems(figure, drawn)
    png = b""
    if not problems:
        png = render_png(figure)
        problems = size_problems(png)
    for problem in problems:
        logger.error(problem)
    if problems:
        return 1
    args.output.write_bytes(png)
    logger.info("wrote %s", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
