# pipeline/scripts/registry.py

import json
from dataclasses import dataclass
from datetime import datetime
from enum import Enum, auto
from pathlib import Path
from typing import Any


class ProcessingStatus(Enum):
    PENDING = auto()
    FAILED = auto()
    DONE = auto()


@dataclass
class StudyFileTracker:
    study_id: str
    file_path: Path
    status: ProcessingStatus = ProcessingStatus.PENDING
    processed_at: datetime | None = None
    error_message: str | None = None

    def complete(self):
        self.status = ProcessingStatus.DONE
        self.processed_at = datetime.now()
        self.error_message = None

    def fail(self, error_message: str):
        self.status = ProcessingStatus.FAILED
        self.error_message = error_message

    def to_dict(self) -> dict[str, Any]:
        return {
            "study_id": self.study_id,
            "file_path": str(self.file_path),
            "status": self.status.name,
            "processed_at": self.processed_at.isoformat()
            if self.processed_at
            else None,
            "error_message": self.error_message,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "StudyFileTracker":
        return cls(
            study_id=data["study_id"],
            file_path=Path(data["file_path"]),
            status=ProcessingStatus[data["status"]],
            processed_at=datetime.fromisoformat(data["processed_at"])
            if data["processed_at"]
            else None,
            error_message=data["error_message"],
        )


class StudyProcessingRegistry:
    def __init__(self):
        self._studies: dict[str, StudyFileTracker] = {}

    def register(self, study_id: str, file_path: Path) -> StudyFileTracker:
        if study_id in self._studies:
            raise ValueError(f"Study '{study_id}' is already registered")
        tracker = StudyFileTracker(study_id=study_id, file_path=file_path)
        self._studies[study_id] = tracker
        return tracker

    def get(self, study_id: str) -> StudyFileTracker:
        try:
            return self._studies[study_id]
        except KeyError:
            raise KeyError(f"Study '{study_id}' not found in registry") from None

    def all_pending(self) -> list[StudyFileTracker]:
        return [
            s for s in self._studies.values() if s.status == ProcessingStatus.PENDING
        ]

    def all_failed(self) -> list[StudyFileTracker]:
        return [
            s for s in self._studies.values() if s.status == ProcessingStatus.FAILED
        ]

    def all_done(self) -> list[StudyFileTracker]:
        return [s for s in self._studies.values() if s.status == ProcessingStatus.DONE]

    @property
    def summary(self) -> dict[str, int]:
        return {
            status.name: sum(1 for s in self._studies.values() if s.status == status)
            for status in ProcessingStatus
        }

    # --- Persistence ---

    def save(self, path: Path) -> None:
        data = [tracker.to_dict() for tracker in self._studies.values()]
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)

    @classmethod
    def load(cls, path: Path) -> "StudyProcessingRegistry":
        """Load the registry from a JSON file. Returns an empty registry if the file does not exist."""
        registry = cls()
        if not path.exists():
            return registry
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        for item in data:
            tracker = StudyFileTracker.from_dict(item)
            registry._studies[tracker.study_id] = tracker
        return registry
