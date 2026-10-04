"""Keep several strategies proposing at once without one of them exiting while the other holds every worker."""
from __future__ import annotations

import asyncio
from typing import Callable, Iterable, List, Set

from search.client import EvalClient


async def run_streaming(client: EvalClient, room: Callable[[Set[asyncio.Future]], bool],
                        fill: Callable[[], List[asyncio.Future]],
                        on_done: Callable[[asyncio.Future], None],
                        initial: Iterable[asyncio.Future] = ()) -> None:
    """Call `fill` whenever `room` says a worker is free. `fill` returns the futures it just proposed
    (empty when it has nothing worth running). Finished futures are handed to `on_done`.

    Yields after every proposal so a second strategy gets the event loop. When this strategy has nothing in
    flight and the executor reports no free worker, it waits for any result (including another strategy's)
    instead of returning and leaving that strategy alone for the rest of the round.
    """
    wake = asyncio.Event()
    client.on_result.append(lambda _outcome: wake.set())
    pending: Set[asyncio.Future] = set(initial)
    while True:
        while not client.exhausted and room(pending):
            new = fill()
            if not new:
                break
            pending.update(new)
            await asyncio.sleep(0)
        if client.exhausted and not pending:
            return
        if not pending:
            idle = client.idle_slots()
            if idle is None or idle > 0 or client.exhausted:
                return
            wake.clear()
            try:
                await asyncio.wait_for(wake.wait(), timeout=2.0)
            except asyncio.TimeoutError:
                continue
            continue
        done, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
        for fut in done:
            on_done(fut)
        if client.exhausted and not pending:
            return
