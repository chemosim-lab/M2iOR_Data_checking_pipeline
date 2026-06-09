import json
import os
import hashlib
from typing import Any

CACHE_DIR = "./cache"


def _path(key: str, subdir: str) -> str:
    # hash pour éviter les caractères interdits dans les noms de fichiers
    safe = hashlib.md5(key.encode()).hexdigest()
    d = os.path.join(CACHE_DIR, subdir)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{safe}.json")


def set_cache(key: str, data: dict[str, Any | None], subdir: str):
    with open(_path(key, subdir), "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def get_cache(key: str, subdir: str) -> Any | None:
    p = _path(key, subdir)
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    return None


def delete_cache(key: str, subdir: str):
    p = _path(key, subdir)
    if os.path.exists(p):
        os.remove(p)


def get_missing_keys(keys: list[str], subdir: str) -> list[str]:
    return [key for key in keys if not os.path.exists(_path(key, subdir))]
