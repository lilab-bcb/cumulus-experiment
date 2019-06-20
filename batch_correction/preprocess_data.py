import scCloud
import sys

def process_data(in_file):
	output_name = in_file.split('.')[0] + "_filter_norm"

	adata = scCloud.tools.read_input(in_file)
	scCloud.tools.update_var_names(adata)
	scCloud.tools.filter_data(adata)
	scCloud.tools.log_norm(adata, norm_count = 1e5)
	scCloud.tools.write_output(adata, output_name)

if __name__ == '__main__':
	in_file = sys.argv[1]
	process_data(in_file)