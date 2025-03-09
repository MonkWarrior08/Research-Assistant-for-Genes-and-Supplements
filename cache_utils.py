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
    
