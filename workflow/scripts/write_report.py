def write_report(report_path, *args):
    with open(report_path, "w") as rep:
        for arg in args:
            rep.write(f"Outfile {arg} was succesfully generated")

if __name__ == "__main__":
    import sys
    report_path = sys.argv[1]
    args = sys.argv[2].strip().split(",")
