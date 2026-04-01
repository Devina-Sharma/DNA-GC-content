from Bio import SeqIO


def gc_content(sequence):
    sequence = sequence.upper()
    g = sequence.count("G")
    c = sequence.count("C")
    return (g + c) / len(sequence) * 100 if len(sequence) > 0 else 0


def process_fasta(file):
    results = []
    for record in SeqIO.parse(file, "fasta"):
        seq = str(record.seq)
        gc = gc_content(seq)
        results.append({
            "id": record.id,
            "length": len(seq),
            "gc_content": gc
        })
    return results


def save_results(results, filename):
    with open(filename, "w") as f:
        f.write("Sequence_ID\tLength\tGC_Content(%)\n")
        for r in results:
            f.write(f"{r['id']}\t{r['length']}\t{r['gc_content']:.2f}\n")


if __name__ == "__main__":
    fasta_file = "sequences.fasta"

    results = process_fasta(fasta_file)

    print("GC CONTENT RESULTS:\n")
    for r in results:
        print(f"{r['id']} | Length: {r['length']} | GC%: {r['gc_content']:.2f}")

    save_results(results, "gc_results.txt")

    print("\nResults saved to gc_results.txt")
