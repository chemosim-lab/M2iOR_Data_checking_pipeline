import json
import os
import hashlib
from typing import Any

CACHE_DIR = "./cache"
os.makedirs(CACHE_DIR, exist_ok=True)


def _path(key: str) -> str:
    # hash pour éviter les caractères interdits dans les noms de fichiers
    safe = hashlib.md5(key.encode()).hexdigest()
    return os.path.join(CACHE_DIR, f"{safe}.json")


def set_cache(key: str, data: dict[str, Any | None]):
    with open(_path(key), "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def get_cache(key: str) -> Any | None:
    p = _path(key)
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    return None


def delete_cache(key: str):
    p = _path(key)
    if os.path.exists(p):
        os.remove(p)


def get_missing_keys(keys: list[str]) -> list[str]:
    return [key for key in keys if not os.path.exists(_path(key))]
