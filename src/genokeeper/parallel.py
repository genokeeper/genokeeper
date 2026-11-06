"""A thread pool's map that stops at the first failure."""

from collections.abc import Callable, Iterable
from concurrent.futures import (
    FIRST_EXCEPTION,
    CancelledError,
    ThreadPoolExecutor,
    wait,
)


def map_until_failure[T, R](
    fn: Callable[[T], R], items: Iterable[T], *, workers: int
) -> list[R]:
    """`fn` over `items`, at most `workers` at a time, the results in `items`' order. At the
    first exception an item no worker has taken never starts, calls already inside `fn` run
    to their end, even one still waiting at the core gate to start its container, and that
    exception is raised. `ThreadPoolExecutor.map` instead goes on through the queue, and
    hides a failure behind every earlier item still running."""
    failures: list[Exception] = []

    def guarded(item: T) -> R:
        # Cancelling is not enough: the worker that failed takes the next item before the
        # waiting thread wakes (200 of 200 trials, CPython 3.14).
        if failures:
            raise CancelledError
        try:
            return fn(item)
        except Exception as e:
            failures.append(e)
            raise

    pool = ThreadPoolExecutor(max_workers=workers)
    try:
        futures = [pool.submit(guarded, item) for item in items]
        wait(futures, return_when=FIRST_EXCEPTION)
    finally:
        pool.shutdown(cancel_futures=True)
    if failures:
        raise failures[0]
    return [future.result() for future in futures]
