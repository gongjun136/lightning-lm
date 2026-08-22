#!/usr/bin/env python3
"""Emit FIRST_POSE/LOST/RECOVERED from one persistent /PosRes subscription."""

from __future__ import annotations

import argparse
import math
import time
from typing import Any

import rclpy
from rclpy.node import Node
from rclpy.qos import DurabilityPolicy, HistoryPolicy, QoSProfile, ReliabilityPolicy


class PosResSilenceMonitor(Node):
    def __init__(self, timeout_sec: float, confirmation_sec: float, message_type: Any) -> None:
        super().__init__("lightning_posres_silence_monitor")
        self.timeout_sec = timeout_sec
        self.confirmation_sec = confirmation_sec
        self.last_arrival: float | None = None
        self.loss_candidate_since: float | None = None
        self.lost = False
        qos = QoSProfile(
            history=HistoryPolicy.KEEP_LAST,
            depth=1,
            reliability=ReliabilityPolicy.BEST_EFFORT,
            durability=DurabilityPolicy.VOLATILE,
        )
        self.create_subscription(message_type, "/PosRes", self.on_posres, qos)
        self.create_timer(min(0.1, timeout_sec / 10.0), self.check_silence)

    @staticmethod
    def emit(event: str, silence_sec: float) -> None:
        print(f"{event},{silence_sec:.6f}", flush=True)

    def on_posres(self, _: Any) -> None:
        now = time.monotonic()
        if self.last_arrival is None:
            self.emit("FIRST_POSE", 0.0)
        elif self.lost:
            self.emit("RECOVERED", max(0.0, now - self.last_arrival))
        self.last_arrival = now
        self.loss_candidate_since = None
        self.lost = False

    def check_silence(self) -> None:
        if self.last_arrival is None or self.lost:
            return
        now = time.monotonic()
        silence_sec = now - self.last_arrival
        if silence_sec < self.timeout_sec:
            self.loss_candidate_since = None
            return
        if self.loss_candidate_since is None:
            # A short confirmation window lets already-queued DDS data run before
            # declaring loss after a temporary CPU scheduling stall.
            self.loss_candidate_since = now
            return
        if now - self.loss_candidate_since >= self.confirmation_sec:
            self.lost = True
            self.emit("LOST", silence_sec)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout-sec", type=float, required=True)
    parser.add_argument("--confirmation-sec", type=float, default=0.25)
    args = parser.parse_args()
    if not math.isfinite(args.timeout_sec) or args.timeout_sec <= 0.0:
        parser.error("--timeout-sec must be finite and positive")
    if not math.isfinite(args.confirmation_sec) or args.confirmation_sec < 0.0:
        parser.error("--confirmation-sec must be finite and non-negative")

    try:
        from geosun_msgs.msg import PosRes
    except ModuleNotFoundError as exc:
        parser.error(f"geosun_msgs is required to monitor /PosRes: {exc}")

    rclpy.init()
    node = PosResSilenceMonitor(args.timeout_sec, args.confirmation_sec, PosRes)
    try:
        rclpy.spin(node)
    except (KeyboardInterrupt, rclpy.executors.ExternalShutdownException):
        pass
    finally:
        node.destroy_node()
        if rclpy.ok():
            rclpy.shutdown()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
