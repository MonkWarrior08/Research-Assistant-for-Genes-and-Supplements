import os
import json
import hashlib
from datetime import datetime, timedelta

cache_dir = ".cache"

def cache_exist():
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

def cache_key(query, func_name):
    key = f"{func_name}_{query}"
    return hashlib.md5(key.encode('utf-8')).hexdigest()

def save_cache(key, data, expiry=7):
    cache_exist()
    cache_file = os.path.join(cache_dir, f"{key}.json")
    expiry_date = (datetime.now() + timedelta(days=expiry)).isoformat()
    cache_data = {
        "data": data,
        "expiry": expiry_date
    }
    with open(cache_file, 'w') as f:
        json.dump(cache_data, f)

def load_cache(key):
    cache_exist()
    cache_file = os.path.join(cache_dir, f"{key}.json")
    if not os.path.exists(cache_file):
        return None
    
    try:
        with open(cache_file, 'r') as f:
            cache_data = json.load(f)
        
        expiry = datetime.fromisoformat(cache_data["expiry"])
        if datetime.now() > expiry:
            os.remove(cache_file)
            return None
        
        return cache_data["data"]

    except (json.JSONDecodeError, KeyError, ValueError):
        if os.path.exists(cache_file):
            os.remove(cache_file)
        return None

def cached_function(func):
    def wrapper(*args, **kwargs):
        if args and isinstance(args[0], str):
            query = args[0]
            key = cache_key(query, func.__name__)
            cache_result = load_cache(key)

            if cache_result is not None:
                return cache_result
            result = func(*args, **kwargs)
            return result
        else:
            return func(*args, **kwargs)
    return wrapper