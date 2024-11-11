c = {'dt': 0.1}

for key, values in c.items():
    exec(f"{(key)} = {values}")

print(dt)


