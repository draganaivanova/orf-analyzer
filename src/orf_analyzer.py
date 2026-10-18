
from Bio import SeqIO
import csv


def find_orfs(sequence):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    for frame in range(3):
        start_positions = []

        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i + 3]

            if codon == start_codon:
                start_positions.append(i)

            elif codon in stop_codons and start_positions:
                end_position = i + 3

                for start_position in start_positions:
                    orf_sequence = sequence[start_position:end_position]
                    length = len(orf_sequence)

                    if length < 9:
                        continue

                    protein = orf_sequence.translate()

                    orfs.append(
                        (
                            frame,
                            orf_sequence,
                            start_position,
                            end_position,
                            length,
                            protein,
                        )
                    )

                start_positions = []

    return orfs


def main():
    all_orfs = []

    for record in SeqIO.parse("data/example_sequences.fasta", "fasta"):
        print("\nSequence ID:", record.id)
        print("Sequence:", record.seq)

        reverse_sequence = record.seq.reverse_complement()
        print("Reverse complement:", reverse_sequence)

        orfs_found = find_orfs(record.seq)
        reverse_orfs = find_orfs(reverse_sequence)
        sequence_length = len(record.seq)

        total_orfs = len(orfs_found) + len(reverse_orfs)
        print("\nTotal ORFs found:", total_orfs)

        # Collect ORFs from the forward strand
        for frame, orf_sequence, start_position, end_position, length, protein in orfs_found:
            print(
                "ORF:", orf_sequence,
                "| Strand: +",
                "| Frame:", frame + 1,
                "| Start:", start_position,
                "| End:", end_position,
                "| Length:", length,
                "| Protein:", protein,
            )

            all_orfs.append({
                "sequence_id": record.id,
                "strand": "+",
                "frame": frame + 1,
                "start": start_position,
                "end": end_position,
                "length": length,
                "orf_sequence": str(orf_sequence),
                "protein": str(protein),
            })

        # Collect ORFs from the reverse strand
        for frame, orf_sequence, start_position, end_position, length, protein in reverse_orfs:
            original_start = sequence_length - end_position
            original_end = sequence_length - start_position

            print(
                "ORF:", orf_sequence,
                "| Strand: -",
                "| Frame:", frame + 1,
                "| Start (original sequence):", original_start,
                "| End (original sequence):", original_end,
                "| Length:", length,
                "| Protein:", protein,
            )

            all_orfs.append({
                "sequence_id": record.id,
                "strand": "-",
                "frame": frame + 1,
                "start": original_start,
                "end": original_end,
                "length": length,
                "orf_sequence": str(orf_sequence),
                "protein": str(protein),
            })

    # Export all detected ORFs to CSV
    with open("results/orfs.csv", "w", newline="") as csvfile:
        fieldnames = [
            "sequence_id",
            "strand",
            "frame",
            "start",
            "end",
            "length",
            "orf_sequence",
            "protein",
        ]

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_orfs)

    print("\nORF results saved to results/orfs.csv")


if __name__ == "__main__":
    main()
