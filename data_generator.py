#!/usr/bin/env python
import serial, time, csv, sys

def main(port='COM4', baud=115200, timeout=5, output_file='datos_planta.csv'):
    try:
        ser = serial.Serial(port, baud, timeout=timeout)
    except Exception as e:
        print(f"Error abriendo {port}: {e}")
        sys.exit(1)

    time.sleep(2)
    ser.write(b'START\n')
    datos = []
    while True:
        line = ser.readline().decode('utf-8', errors='ignore').strip()
        if not line: continue
        if line == 'END': break
        if line.startswith('t,'): continue
        parts = line.split(',')
        if len(parts)!=4: continue
        try:
            t, u, v, raw = float(parts[0]), float(parts[1]), float(parts[2]), int(parts[3])
        except ValueError:
            continue
        datos.append([t,u,v,raw])
        print(f"{t:.3f}, {u:.3f}, {v:.3f}, {raw}")

    ser.close()
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['t','u','v','raw'])
        writer.writerows(datos)
    print(f"Guardado {len(datos)} muestras en '{output_file}'")

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('--port',default='COM4')
    p.add_argument('--baud',type=int,default=115200)
    p.add_argument('--output',default='datos_planta.csv')
    args=p.parse_args()
    main(args.port,args.baud,5,args.output)
