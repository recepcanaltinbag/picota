"""
Shared fixtures for picota unit tests.
"""
import sys
import os

# picota/picota/ dizinini path'e ekle — kaynak kod "from src.X import" kullanıyor
tests_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(tests_dir)            # picota/
inner_root   = os.path.join(project_root, 'picota')  # picota/picota/
sys.path.insert(0, inner_root)
