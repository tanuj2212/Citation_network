import requests
import time
import pandas as pd
import xml.etree.ElementTree as ET
import pickle  # To save citation network
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

# Timer start
start_time = time.time()

# PubMed API URLs
PUBMED_SEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_SUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
PUBMED_CITATION_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"

# COVID-related MeSH terms
MESH_TERMS = [
    "COVID-19", "SARS-CoV-2 Infection", "Coronavirus Infections",
    "Severe Acute Respiratory Syndrome (SARS)", "Post-Acute COVID-19 Syndrome (Long COVID)",
    "COVID-19 Vaccines", "COVID-19 Testing (PCR, Antigen)", "Lockdowns & Social Distancing",
    "PPE (Personal Protective Equipment)", "Novel Coronavirus", "Wuhan Coronavirus",
    "Coronavirus Disease 2019"
]

# Constants
CURRENT_YEAR = 2025
LIFETIME = 5  # Years


# Function to fetch all COVID-19 papers from PubMed (2019-2024)
def fetch_pubmed_ids():
    query = " OR ".join([f"{term}[MeSH]" for term in MESH_TERMS])
    params = {
        "db": "pubmed",
        "term": f"({query}) AND 2019:2024[PDAT]",
        "retmode": "xml",
        "retmax": 500000
    }
    response = requests.get(PUBMED_SEARCH_URL, params=params)
    root = ET.fromstring(response.content)
    return [id_elem.text for id_elem in root.findall(".//Id")]


# Function to fetch metadata (Title, DOI, Publication Date)
def fetch_pubmed_metadata(pmids):
    metadata = {}
    batch_size = 50
    max_retries = 5

    def fetch_batch(batch_pmids):
        params = {"db": "pubmed", "id": ",".join(batch_pmids), "retmode": "json"}

        for attempt in range(max_retries):
            response = requests.get(PUBMED_SUMMARY_URL, params=params)
            if response.status_code == 429:
                time.sleep(2 ** attempt)
                continue
            if response.status_code == 200:
                try:
                    results = response.json()["result"]
                    return {pmid: {
                        "Title": results[pmid].get("title", "N/A"),
                        "Publication Date": results[pmid].get("pubdate", "N/A"),
                        "DOI": results[pmid].get("elocationid", "N/A")
                    } for pmid in batch_pmids if pmid in results}
                except:
                    pass
            break
        return {}

    with ThreadPoolExecutor(max_workers=5) as executor:
        batches = [pmids[i:i + batch_size] for i in range(0, len(pmids), batch_size)]
        results = executor.map(fetch_batch, batches)

    for result in results:
        metadata.update(result)

    return metadata


# Function to fetch citations & references
import requests
import time
import json
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from requests.exceptions import ChunkedEncodingError, RequestException

PUBMED_CITATION_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"

def fetch_citations_and_references(pmids):
    """Fetch citations and references with error handling and retries."""
    citation_network = defaultdict(list)
    batch_size = 20  # Reduce batch size to avoid errors
    max_retries = 5
    failed_pmids = []  # Store failed PMIDs

    def fetch_batch(batch_pmids):
        """Fetch references and citations for a batch of PMIDs with retries."""
        local_network = defaultdict(list)
        params = {
            "dbfrom": "pubmed", "db": "pubmed",
            "id": ",".join(batch_pmids), "cmd": "neighbor",
            "linkname": "pubmed_pubmed_refs", "retmode": "json"
        }

        for attempt in range(max_retries):
            try:
                response = requests.get(PUBMED_CITATION_URL, params=params, timeout=10, stream=True)
                if response.status_code == 429:
                    wait_time = 2 ** attempt
                    print(f"⚠️ Too many requests! Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                    continue  # Retry request

                if response.status_code == 200:
                    refs_data = response.json()
                    for linkset in refs_data.get("linksets", []):
                        pmid = linkset.get("ids", ["N/A"])[0]
                        if "linksetdbs" in linkset:
                            for linksetdb in linkset["linksetdbs"]:
                                if linksetdb["linkname"] == "pubmed_pubmed_refs":
                                    local_network[pmid].extend(linksetdb.get("links", []))
                    return local_network  # Return successful response

            except (ChunkedEncodingError, RequestException, json.JSONDecodeError) as e:
                print(f"⚠️ Error fetching batch: {e}")
                time.sleep(2)  # Small delay before retrying

        failed_pmids.extend(batch_pmids)  # Store failed PMIDs
        return {}

    # Use parallel processing to fetch citation data
    with ThreadPoolExecutor(max_workers=5) as executor:
        batches = [pmids[i:i + batch_size] for i in range(0, len(pmids), batch_size)]
        results = executor.map(fetch_batch, batches)

    # Merge all results
    for result in results:
        for key, value in result.items():
            citation_network[key].extend(value)

    # Log failed PMIDs
    if failed_pmids:
        with open("failed_pmids.txt", "w") as f:
            f.write("\n".join(failed_pmids))
        print(f"⚠️ {len(failed_pmids)} PMIDs failed. Saved to 'failed_pmids.txt' for retry.")

    return citation_network


# Main Execution
print("\n🔍 Fetching PubMed Papers...")
pmids = fetch_pubmed_ids()
print(f"✅ Retrieved {len(pmids)} COVID-19 papers.")

print("\n📄 Fetching Metadata for Papers...")
metadata = fetch_pubmed_metadata(pmids)

print("\n🔗 Fetching Citation Network...")
citation_network = fetch_citations_and_references(pmids)

# Save citation network as a dump file
with open("covid_citation_network.pkl", "wb") as f:
    pickle.dump(citation_network, f)
print("✅ Citation network saved as dump file.")

# Compute Degree of Connectivity (DC) and Freshness Factor (FF)
dc = {pmid: len(citation_network.get(pmid, [])) for pmid in pmids}
ff = {}
missing_pmids = []

for pmid in pmids:
    if pmid in metadata and metadata[pmid]["Publication Date"] != "N/A":
        ff[pmid] = max(0, (int(metadata[pmid]["Publication Date"][:4]) - 2019) / (datetime.now().year - 2019))
    else:
        ff[pmid] = 0  # Assign default FF value if missing
        missing_pmids.append(pmid)

# Log missing PMIDs for later debugging
if missing_pmids:
    with open("missing_pmids.txt", "w") as f:
        f.write("\n".join(missing_pmids))
    print(f"⚠️ {len(missing_pmids)} PMIDs missing metadata. Logged in 'missing_pmids.txt'.")
# Compute DCIF Score
dcif = {pmid: ff[pmid] * dc[pmid] for pmid in pmids}

# Prepare ranking data
ranking_data = [{
    "PMID": pmid,
    "Title": metadata[pmid]["Title"],
    "Publication Date": metadata[pmid]["Publication Date"],
    "DOI": metadata[pmid]["DOI"],
    "Degree of Connectivity (DC)": dc[pmid],
    "Freshness Factor (FF)": ff[pmid],
    "DCIF Score": dcif[pmid]
} for pmid in pmids]

# Save ranking results
ranking_df = pd.DataFrame(ranking_data)
ranking_df.sort_values(by="DCIF Score", ascending=False, inplace=True)

# Save in CSV
ranking_df.to_csv("covid_ranking.csv", index=False)
print("\n✅ Ranking saved to covid_ranking.csv")

# Timer end
end_time = time.time()
print(f"\n⏳ Total Execution Time: {end_time - start_time:.2f} seconds.")
