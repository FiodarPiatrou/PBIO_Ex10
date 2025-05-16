from Bio import Entrez
import time
import os


class NCBIRetriever:
    def __init__(self, email, api_key):
        """Initialize with NCBI credentials."""
        self.email = email
        self.api_key = api_key

        # Ustawienia Entrez
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid,min_len=None,max_len=None):
        """Search for all records associated with a taxonomic ID."""
        print(f"Searching for records with taxID: {taxid}")
        try:
            is_min_len=min_len is not None
            is_max_len=max_len is not None
            if is_max_len and is_min_len and min_len > max_len:
                print("min_len must be < max_len")
                return None
            # Najpierw pobierz informacje taksonomiczne
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            # Szukaj rekordów
            search_term = f"txid{taxid}[Organism]"
            if is_min_len or is_max_len:
                search_term += f" AND {min_len if is_min_len else 0}:{max_len if is_max_len else '*'}[SLEN]"

            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print(f"No records found for {organism_name}")
                return None

            print(f"Found {count} records")

            # Zapisz wyniki wyszukiwania do późniejszego wykorzystania
            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count

            return count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=10):
        """Fetch a batch of records using the stored search results."""
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            # Limit, aby zapobiec przeciążeniu serwera
            batch_size = min(max_records, 500)

            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )

            # Surowy rekord GenBank
            records_text = handle.read()

            return records_text

        except Exception as e:
            print(f"Error fetching records: {e}")
            return ""


def main():
    # Uzyskaj dane uwierzytelniające
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")

    # Utwórz obiekt retriever
    retriever = NCBIRetriever(email, api_key)

    # Uzyskaj taxid od użytkownika
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")
    min_len = input("Enter minimum length of sequences: ")
    max_len = input("Enter maximum length of sequences: ")
    min_len = int(min_len) if min_len else None
    max_len = int(max_len) if max_len else None
    # Szukaj rekordów
    count = retriever.search_taxid(taxid, min_len, max_len)

    if not count:
        print("No records found. Exiting.")
        return

    # Pobierz kilka pierwszych rekordów jako próbkę
    print("\nFetching sample records...")
    sample_records = retriever.fetch_records(start=0, max_records=5)

    # Zapisz do pliku
    output_file = f"taxid_{taxid}_sample.gb"
    with open(output_file, "w") as f:
        f.write(sample_records)

    print(f"Saved sample records to {output_file}")
    print("\nNote: This is just a basic retriever. You need to extend its functionality!")


if __name__ == "__main__":
    main()