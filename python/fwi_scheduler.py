#!/usr/bin/env python3
"""
Seismic FWI Task Scheduler - Python orchestration layer for the FWI engine.
Provides high-level workflow management for seismic forward modeling tasks.
"""

import os
import sys
import time
import logging
from typing import List, Dict, Optional, Callable, Any
from dataclasses import dataclass, field
from pathlib import Path

try:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'python'))
    import pyfwi as fwi
except ImportError:
    print("Warning: pyfwi module not found. Please build the project first.")
    fwi = None

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class TaskConfig:
    """Configuration for a forward modeling task."""
    nx: int = 200
    nz: int = 150
    dx: float = 10.0
    dz: float = 10.0
    num_sources: int = 5
    num_receivers: int = 100
    peak_frequency: float = 30.0
    num_time_steps: int = 1000
    dt: float = 0.0
    pml_width: int = 20
    snapshot_interval: int = 50
    water_velocity: float = 1500.0
    sediment_velocity: float = 2500.0
    basement_velocity: float = 3500.0
    output_dir: str = "./output"
    vtk_format: str = "ASCII"
    compute_device: str = "CPU"
    shot_indices: List[int] = field(default_factory=lambda: [0])


class TaskStatus:
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class FWITask:
    """Represents a single FWI modeling task."""

    def __init__(self, config: TaskConfig, task_id: str = "task_0"):
        self.config = config
        self.task_id = task_id
        self.status = TaskStatus.PENDING
        self.progress = 0.0
        self.start_time: Optional[float] = None
        self.end_time: Optional[float] = None
        self.error: Optional[str] = None
        self.result: Dict[str, Any] = {}

    def run(self, engine: Optional[fwi.FwiEngine] = None) -> bool:
        if fwi is None:
            raise RuntimeError("pyfwi module not available. Build the C++ library first.")

        self.status = TaskStatus.RUNNING
        self.start_time = time.time()
        logger.info(f"Starting task {self.task_id}...")

        try:
            if engine is None:
                engine = fwi.FwiEngine()

            device = fwi.ComputeDevice.CPU if self.config.compute_device == "CPU" else fwi.ComputeDevice.CUDA
            engine.set_compute_device(device)
            logger.info(f"Using compute device: {self.config.compute_device}")

            engine.create_synthetic_velocity_model(
                self.config.nx, self.config.nz,
                self.config.dx, self.config.dz,
                self.config.water_velocity,
                self.config.sediment_velocity,
                self.config.basement_velocity
            )
            logger.info(f"Created velocity model: {self.config.nx}x{self.config.nz} "
                       f"grid, dx={self.config.dx}m, dz={self.config.dz}m")

            domain_width = (self.config.nx - 1) * self.config.dx
            pml_extent = self.config.pml_width * self.config.dx
            source_start = pml_extent + 50
            source_end = domain_width - pml_extent - 50
            receiver_start = pml_extent + 20
            receiver_end = domain_width - pml_extent - 20

            engine.create_survey(
                source_start, source_end, self.config.num_sources,
                receiver_start, receiver_end, self.config.num_receivers,
                depth=2 * self.config.dz
            )
            logger.info(f"Created survey: {self.config.num_sources} sources, "
                       f"{self.config.num_receivers} receivers")

            wavelet = fwi.Wavelet.ricker(0.001, self.config.peak_frequency)
            engine.set_wavelet(wavelet)
            logger.info(f"Set wavelet: Ricker, {self.config.peak_frequency}Hz peak frequency")

            params = fwi.FDSimulationParams()
            params.num_time_steps = self.config.num_time_steps
            params.dt = self.config.dt
            params.pml_width = self.config.pml_width
            params.snapshot_interval = self.config.snapshot_interval
            params.use_pml = True
            params.record_receivers = True
            engine.set_simulation_params(params)

            vtk_format_map = {
                "ASCII": fwi.VtkFormat.LEGACY_ASCII,
                "BINARY": fwi.VtkFormat.LEGACY_BINARY,
                "XML": fwi.VtkFormat.XML_IMAGE_DATA
            }
            engine.set_vtk_format(vtk_format_map.get(self.config.vtk_format,
                                                    fwi.VtkFormat.LEGACY_ASCII))

            output_dir = Path(self.config.output_dir) / self.task_id
            output_dir.mkdir(parents=True, exist_ok=True)

            vel_file = output_dir / "velocity_model.vtk"
            engine.export_velocity_model(str(vel_file))
            logger.info(f"Exported velocity model to {vel_file}")

            for shot_idx in self.config.shot_indices:
                logger.info(f"Processing shot {shot_idx}...")
                success = engine.run_forward_modeling(shot_idx)
                if not success:
                    raise RuntimeError(f"Forward modeling failed for shot {shot_idx}")

                shot_dir = output_dir / f"shot_{shot_idx}"
                shot_dir.mkdir(exist_ok=True)

                snap_dir = shot_dir / "snapshots"
                snap_dir.mkdir(exist_ok=True)
                engine.export_all_snapshots(str(snap_dir), "snapshot")
                logger.info(f"Exported {len(engine.snapshots())} snapshots to {snap_dir}")

                rec_file = shot_dir / "receiver_data.csv"
                engine.export_receiver_data(str(rec_file))
                logger.info(f"Exported receiver data to {rec_file}")

                if len(engine.snapshots()) > 0:
                    for i, snap in enumerate(engine.snapshots()):
                        if i % 5 == 0:
                            vtk_3d = shot_dir / f"snapshot_3d_{i:06d}.vtk"
                            writer = fwi.VtkWriter()
                            writer.write_3d_snapshot(str(vtk_3d), snap, engine.velocity_model().grid())

            self.status = TaskStatus.COMPLETED
            self.end_time = time.time()
            self.result["output_dir"] = str(output_dir)
            self.result["num_snapshots"] = len(engine.snapshots())
            self.result["num_receivers"] = len(engine.receiver_data())
            self.progress = 1.0

            logger.info(f"Task {self.task_id} completed successfully in "
                       f"{self.end_time - self.start_time:.2f}s")
            return True

        except Exception as e:
            self.status = TaskStatus.FAILED
            self.error = str(e)
            self.end_time = time.time()
            logger.error(f"Task {self.task_id} failed: {e}")
            return False


class FWIScheduler:
    """Manages and schedules multiple FWI tasks."""

    def __init__(self, max_parallel: int = 1):
        self.max_parallel = max_parallel
        self.tasks: List[FWITask] = []
        self.completed_tasks: List[FWITask] = []
        self.failed_tasks: List[FWITask] = []

    def add_task(self, config: TaskConfig) -> str:
        task_id = f"task_{len(self.tasks)}"
        task = FWITask(config, task_id)
        self.tasks.append(task)
        logger.info(f"Added task {task_id}")
        return task_id

    def run_all(self) -> bool:
        logger.info(f"Starting scheduler with {len(self.tasks)} tasks...")
        all_success = True

        for i, task in enumerate(self.tasks):
            logger.info(f"Running task {i+1}/{len(self.tasks)}: {task.task_id}")
            success = task.run()
            if success:
                self.completed_tasks.append(task)
            else:
                self.failed_tasks.append(task)
                all_success = False

        logger.info(f"Scheduler finished. "
                   f"Completed: {len(self.completed_tasks)}, "
                   f"Failed: {len(self.failed_tasks)}")
        return all_success

    def get_status(self) -> Dict[str, int]:
        return {
            "total": len(self.tasks),
            "completed": len(self.completed_tasks),
            "failed": len(self.failed_tasks),
            "pending": len(self.tasks) - len(self.completed_tasks) - len(self.failed_tasks)
        }


def run_forward_modeling_example():
    """Run a complete forward modeling example."""
    if fwi is None:
        print("Error: pyfwi module not available. Build the project first.")
        return False

    print("=" * 60)
    print("Seismic FWI Engine - Forward Modeling Example")
    print("=" * 60)

    config = TaskConfig(
        nx=150,
        nz=120,
        dx=10.0,
        dz=10.0,
        num_sources=3,
        num_receivers=80,
        peak_frequency=25.0,
        num_time_steps=800,
        dt=0.0008,
        pml_width=15,
        snapshot_interval=40,
        output_dir="./output",
        vtk_format="ASCII",
        compute_device="CPU",
        shot_indices=[0, 1]
    )

    scheduler = FWIScheduler()
    scheduler.add_task(config)
    success = scheduler.run_all()

    if success:
        print("\nForward modeling completed successfully!")
        print(f"Output directory: {config.output_dir}")
    else:
        print("\nForward modeling failed!")

    return success


def analyze_results(output_dir: str = "./output/task_0"):
    """Analyze and visualize results from a completed simulation."""
    import numpy as np

    output_path = Path(output_dir)
    if not output_path.exists():
        print(f"Output directory {output_dir} does not exist!")
        return

    print("\n" + "=" * 60)
    print("Simulation Results Analysis")
    print("=" * 60)

    if fwi is not None:
        print("\nAvailable wavelet types:")
        for wt in [fwi.WaveletType.RICKER, fwi.WaveletType.GAUSSIAN,
                  fwi.WaveletType.ORMSBY, fwi.WaveletType.SINE]:
            print(f"  - {wt}")

        print("\nAvailable VTK formats:")
        for vf in [fwi.VtkFormat.LEGACY_ASCII, fwi.VtkFormat.LEGACY_BINARY,
                  fwi.VtkFormat.XML_IMAGE_DATA]:
            print(f"  - {vf}")

    print(f"\nOutput contents:")
    for item in sorted(output_path.rglob("*")):
        if item.is_file():
            rel_path = item.relative_to(output_path)
            size_mb = item.stat().st_size / (1024 * 1024)
            print(f"  {rel_path}: {size_mb:.2f} MB")

    print("\nTo view the results:")
    print("  1. Open ParaView or Visit")
    print("  2. Load the velocity_model.vtk file")
    print("  3. Load the snapshot series from snapshots/ directory")
    print("  4. Apply the 'Warp By Scalar' filter for better visualization")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Seismic FWI Engine Task Scheduler")
    parser.add_argument("--run", action="store_true", help="Run forward modeling example")
    parser.add_argument("--analyze", type=str, help="Analyze results from output directory")
    parser.add_argument("--config", type=str, help="Path to JSON config file")

    args = parser.parse_args()

    if args.run:
        run_forward_modeling_example()
    elif args.analyze:
        analyze_results(args.analyze)
    else:
        parser.print_help()
