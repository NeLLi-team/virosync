"""Tests for the README workflow-map generator and the checks it runs before writing."""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import pytest
from PIL import Image

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import plot_workflow_map as wm  # noqa: E402


def _needs_font() -> None:
    try:
        wm.font_family()
    except FileNotFoundError:
        pytest.skip(f"{wm.FONT} is not installed")


def _stop(label_start: str) -> wm.Stop:
    return next(stop for stop in wm.STOPS if stop.label.startswith(label_start))


def _moved(label_start: str, **changes: object) -> tuple[wm.Stop, ...]:
    return tuple(
        dataclasses.replace(stop, **changes) if stop.label.startswith(label_start) else stop for stop in wm.STOPS
    )


def _shifted(track: wm.Track, y: float) -> tuple[wm.Route, ...]:
    return tuple(
        dataclasses.replace(route, points=tuple((x, y) for x, _ in route.points)) if route.track == track else route
        for route in wm.ROUTES
    )


def _annotation_branch_from(y: float) -> tuple[wm.Route, ...]:
    def replace(route: wm.Route) -> wm.Route:
        if route.track != wm.ANNOTATION:
            return route
        return dataclasses.replace(route, points=((wm.ANNOTATION_RISER_X, y), *route.points[2:]))

    return tuple(replace(route) for route in wm.ROUTES)


def test_topology_problems_accepts_the_shipped_map() -> None:
    assert wm.topology_problems(wm.STOPS, wm.ROUTES) == []


def test_topology_problems_compares_lines_by_value() -> None:
    routes = tuple(dataclasses.replace(route, track=dataclasses.replace(route.track)) for route in wm.ROUTES)

    assert wm.topology_problems(wm.STOPS, routes) == []


def test_topology_problems_reports_a_stop_off_a_line_it_lists() -> None:
    stops = _moved("Host\nsignature\nmodel", y_low=wm.Y_HOST - 20, y_high=wm.Y_HOST - 20)

    problems = wm.topology_problems(stops, wm.ROUTES)

    assert any("misses the host line" in problem for problem in problems)


def test_topology_problems_reports_a_line_through_a_stop_that_does_not_list_it() -> None:
    stops = _moved("Composition", y_low=wm.Y_HOST, kind="interchange")

    problems = wm.topology_problems(stops, wm.ROUTES)

    assert any("host line runs through 'Composition" in problem for problem in problems)


def test_topology_problems_reports_a_line_that_only_passes_a_stop_corner() -> None:
    stop = wm.Stop("probe", 180.0, 170.0, 170.0, (wm.ANNOTATION,), "above")
    route = wm.Route(wm.ANNOTATION, ((185.5, 175.5), (195.0, 185.0)), None, None)

    problems = wm.topology_problems([stop], [route])

    assert "'probe' misses the annotation line" in problems


def test_topology_problems_reports_stops_that_touch() -> None:
    stops = _moved("Composition", x=_stop("Gap merge").x + 6)

    problems = wm.topology_problems(stops, wm.ROUTES)

    assert any(problem.startswith("stops ") for problem in problems)


def test_topology_problems_reports_lines_that_touch() -> None:
    routes = _shifted(wm.HOST, wm.Y_MAIN - 4)

    problems = wm.topology_problems(wm.STOPS, routes)

    assert any("host line and the main line touch" in problem for problem in problems)


def test_topology_problems_reports_a_branch_detached_from_its_line() -> None:
    routes = _annotation_branch_from(wm.Y_MAIN + 17)

    problems = wm.topology_problems(wm.STOPS, routes)

    assert "the annotation line is not attached to the main line" in problems


def test_topology_problems_reports_a_branch_from_the_wrong_line() -> None:
    routes = _annotation_branch_from(wm.Y_HOST)

    problems = wm.topology_problems(wm.STOPS, routes)

    assert {"the host line and the annotation line touch", "the annotation line and the host line touch"} & set(
        problems
    )


def test_topology_problems_reports_a_riser_on_a_band_edge() -> None:
    edge = wm.PHASES[2].x1
    routes = tuple(
        dataclasses.replace(
            route, points=tuple((edge, y) if x == wm.ANNOTATION_RISER_X else (x, y) for x, y in route.points)
        )
        if route.track == wm.ANNOTATION
        else route
        for route in wm.ROUTES
    )

    problems = wm.topology_problems(wm.STOPS, routes)

    assert any("runs along the edge of phase" in problem for problem in problems)


def test_riser_problems_ignore_a_line_below_the_bands() -> None:
    route = wm.Route(wm.MAIN, ((wm.PHASES[1].x0, 4.0), (wm.PHASES[1].x0, 18.0)), None, None)

    assert wm.riser_problems([route]) == []


def test_layout_problems_accepts_the_shipped_map() -> None:
    _needs_font()

    figure, drawn = wm.draw_map()

    assert wm.layout_problems(figure, drawn) == []


@pytest.mark.parametrize(
    ("label_start", "changes", "expected"),
    [
        ("Composition", {"x": 400.0}, "too close"),
        ("Host\nsignature\ntrim", {"x": 266.0}, "leaves phase 2"),
        ("Domains", {"x": 535.0}, "nearer to"),
    ],
)
def test_layout_problems_reports_a_misplaced_label(label_start: str, changes: dict, expected: str) -> None:
    _needs_font()

    figure, drawn = wm.draw_map(_moved(label_start, **changes))

    assert any(expected in problem for problem in wm.layout_problems(figure, drawn))


def test_layout_problems_reports_text_below_the_readable_size(monkeypatch: pytest.MonkeyPatch) -> None:
    _needs_font()
    monkeypatch.setattr(wm, "MIN_TEXT_PX", 13.0)

    figure, drawn = wm.draw_map()

    assert any("renders at" in problem for problem in wm.layout_problems(figure, drawn))


def test_main_writes_a_png_of_the_planned_size(tmp_path: Path) -> None:
    _needs_font()
    output = tmp_path / "map.png"

    assert wm.main(["--output", str(output)]) == 0

    with Image.open(output) as image:
        assert image.size == wm.planned_size()


def test_main_keeps_the_old_png_when_the_size_check_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _needs_font()
    output = tmp_path / "map.png"
    output.write_bytes(b"old map")
    monkeypatch.setattr(wm, "planned_size", lambda: (1, 1))

    assert wm.main(["--output", str(output)]) == 1

    assert output.read_bytes() == b"old map"
