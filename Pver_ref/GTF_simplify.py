import re
import argparse

def process_gtf(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            attributes = fields[8]
            
            # Extract required attributes
            gene_id_match = re.search(r'gene_id "[^"]+"', attributes)
            transcript_id_match = re.search(r'transcript_id "[^"]+"', attributes)
            gene_biotype_match = re.search(r'gene_biotype "[^"]+"', attributes)
            
            # Build new attributes string
            new_attributes = []
            if gene_id_match:
                new_attributes.append(gene_id_match.group())
            if transcript_id_match:
                new_attributes.append(transcript_id_match.group())
            if gene_biotype_match:
                new_attributes.append(gene_biotype_match.group())
            
            new_attributes_str = '; '.join(new_attributes)
            if not new_attributes_str.endswith(';'):
                new_attributes_str += ';'
            
            fields[8] = new_attributes_str
            outfile.write('\t'.join(fields) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simplify GTF file attributes.")
    parser.add_argument("input_file", help="Path to the input GTF file")
    parser.add_argument("output_file", help="Path to the output GTF file")
    args = parser.parse_args()
    
    process_gtf(args.input_file, args.output_file)
