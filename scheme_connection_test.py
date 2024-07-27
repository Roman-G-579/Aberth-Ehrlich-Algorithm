import subprocess

result = subprocess.run(['racket', 'utils.rkt', str(3), str(5)], capture_output=True, text=True)
print(result.stdout)


# def call_racket_add(a, b):
#     result = subprocess.run(['racket', 'script.rkt', str(a), str(b)], capture_output=True, text=True)
#     return result.stdout.strip()
#
#
# print(call_racket_add(3, 5))