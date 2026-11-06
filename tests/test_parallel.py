"""`map_until_failure`: the order of its results, and what a failure stops."""

import threading

import pytest

from genokeeper.parallel import map_until_failure


def test_results_keep_the_input_order() -> None:
    second_done = threading.Event()

    def fn(item: int) -> int:
        if item == 0:
            assert second_done.wait(timeout=5)
        else:
            second_done.set()
        return item * 10

    assert map_until_failure(fn, [0, 1], workers=2) == [0, 10]


def test_a_later_failure_stops_what_has_not_started() -> None:
    """Item 1 fails while item 0 still runs: items 2 on never start, item 0 finishes, and
    the failure reaches the caller. On `list(pool.map(...))` the workers go on to item 2,
    which ends item 0's wait, and to every item after it."""
    started: list[int] = []
    finished: list[int] = []
    first_running = threading.Event()
    queued_item_ran = threading.Event()

    def fn(item: int) -> int:
        started.append(item)
        if item == 0:
            first_running.set()
            queued_item_ran.wait(timeout=0.2)
            finished.append(item)
        elif item == 1:
            assert first_running.wait(timeout=5)
            raise ValueError("item 1 failed")
        else:
            queued_item_ran.set()
        return item

    with pytest.raises(ValueError, match="item 1 failed"):
        map_until_failure(fn, range(6), workers=2)
    assert sorted(started) == [0, 1]
    assert finished == [0]
