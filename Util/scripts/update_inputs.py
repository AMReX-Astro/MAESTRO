# a script to replace old-style BC parameters with the new
# string-based ones

import os

# old : new parameter names
bc_param_names = {
    "bcx_lo": "xlo_boundary_type",
    "bcx_hi": "xhi_boundary_type",
    "bcy_lo": "ylo_boundary_type",
    "bcy_hi": "yhi_boundary_type",
    "bcz_lo": "zlo_boundary_type",
    "bcz_hi": "zhi_boundary_type"}

bc_value_names = {
    "11": "inlet",
    "12": "outlet",
    "13": "symmetry",
    "14": "slip wall",
    "15": "no slip wall",
    "-1": "periodic"}



def update_bcs(files):

    for ifile in files:

        new_lines = []
        update = False

        # check to see if any of our BC parameters are in the file,
        # if so, update their format
        with open(ifile, "r") as fin:
            for line in fin:
                nline = None
                for k in bc_param_names:
                    if k in line:
                        update = True
                        ok, ov = line.split("=")
                        if ov.strip() in bc_value_names:
                            nk = bc_param_names[k]
                            nv = '"{}"'.format(bc_value_names[ov.strip()])
                            nline = " = ".join([nk, nv]) + "\n"
                if nline is None:
                    nline = line
                new_lines.append(nline)
                        
        # write out a temporary file, starting with _inputs...
        if update:
            dir_name = os.path.dirname(ifile)
            fname = "_{}".format(os.path.basename(ifile))
            ofile = "/".join([dir_name, fname])

            with open(ofile, "w") as fout:
                for line in new_lines:
                    fout.write(line)

            # mv the new file to the old file
            os.rename(ofile, ifile)


def find_inputs():
    """ find all of the inputs files under the current directory """

    top_dir = os.getcwd()

    inputs_files = []

    for root, sub_folders, files in os.walk(top_dir):
        for f in files:
            if f.startswith("inputs"):
                inputs_files.append(os.path.join(root, f))

    return inputs_files


if __name__ == "__main__":
    inputs = find_inputs()

    update_bcs(inputs)


