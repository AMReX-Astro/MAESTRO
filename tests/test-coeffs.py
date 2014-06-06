
import subprocess

# XXX: these tests don't really do anything yet...

def test_coeffs():
    subprocess.check_call(["tests/coeffs.exe"])

def test_t1():
    subprocess.check_call(["tests/t1.exe"])
