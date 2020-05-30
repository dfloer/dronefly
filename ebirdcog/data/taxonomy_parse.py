import csv
import re
import json
from collections import defaultdict


class Taxonomy:
    """
    Stores a single species' information.
    """

    def __init__(self, common_name, scientific_name, species_code, short_codes):
        self.common_name = common_name
        self.scientific_name = scientific_name
        self.species_code = species_code
        self.short_codes = short_codes

    def __str__(self):
        return f"{self.common_name}, {self.scientific_name}, {self.species_code}, {self.short_codes}"

    def __repr__(self):
        return str(
            {
                "common_name": self.common_name,
                "scientific_name": self.scientific_name,
                "species_code": self.species_code,
                "short_codes": self.short_codes,
            }
        )

    def as_dict(self):
        return {
                "common_name": self.common_name,
                "scientific_name": self.scientific_name,
                "species_code": self.species_code,
                "short_codes": self.short_codes,
            }


def open_raw_csv_ebird(csv_path):
    """
    Opens a raw csv file and parses it into a list of CSV lines.
    Args:
        csv_path (str): Path to the file to open.
    Returns:
        A list of lists, where each sublist is a line from the opened CSV file.
    """
    output = []
    with open(csv_path, "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        raw_lines = list(csv_reader)
    header = raw_lines[0]
    keys = [x.lower() for x in header]
    for entry in raw_lines:
        cells = {}
        for idx, field in enumerate(entry):
            cells[keys[idx]] = field
        output += [cells]
    return output


def name_to_4lc(name, alt_hyphen_split=False):
    """
    Converts a given name to a 4-letter code.
    Uses rules at: https://support.ebird.org/en/support/solutions/articles/48000960508-ebird-mobile-tips-tricks
    The rules don't specifically say how to deal with names shorter than 4 letters, so the whole name is returned in this case.
    Note that eBird's own website does not support searching for some of the birds their own rules generate.
    The are also several ambiguities and edge cases in eBrid's rules for finding 4-letter-codes. A best guess at correct behaviour was taken.
    eBird also used to claim support for both splitting and not-splitting on hyphens in the rules, though this never appeared to be implemented on the website or app.
    Args:
        name (str): name to convert to a 4-letter code
    Returns:
        A list of strings of 4 letters representing the ebird short code. Includes the base code as well as eBird's alternatives.
    """
    res = set()
    skipped_words = ["of", "and", "the"]
    split_name = tuple(re.split(r"[ -]", name))
    res.add(words_to_code(split_name))
    # Find the alternative for names longer than 4 words.
    if len(split_name) > 4:
        new_name = tuple(e for e in split_name if e not in skipped_words)
        res.add(words_to_code(new_name))
    return list(res)


def words_to_code(split_name):
    """
    Takes a tuple of the words that make up a bird's name and returns the actual 4-letter code using eBird's rules.
    Args:
        split_name (tuple): Split words of a bird's name.
    Returns:
        A string of at most 4 characters containing the 4-letter-code for the given split name.
    """
    res = ""
    if len(split_name) == 1:
        res = split_name[0][0:4].upper()
    elif len(split_name) == 2:
        res = (split_name[0][0:2] + split_name[1][0:2]).upper()
    elif len(split_name) == 3:
        res = (split_name[0][0:1] + split_name[1][0:1] + split_name[2][0:2]).upper()
    elif len(split_name) >= 4:
        res = "".join([split_name[x][0] for x in range(4)]).upper()
    return res


def taxonomy_parse(csv_path):
    """
    Parses the taxonomy csv into four-letter codes (both scientific, banding and common name),
        as well as eBird's unique codes and scientific + common names.
    Takes the raw eBird axonomy csv (not Clements or combined eBird/Clements).
    Note that there are collisions in 4 letter codes. This code does not attempt to disambiguate them.
    Args:
        csv_path (str): Path to the file to open.
    Returns:
        4 dictionaries:
        1. key: each 4 letter code.
        2. key: ebird unique codes.
        3. key: common name.
        4. key:scientific name.
        With the values being an object containing the values associated with that key.
        Why do it this way? To make it easy to look things up by any of the 4 possible key types.
    """
    common_map = {}
    scientific_map = {}
    code_map = {}
    short_map = defaultdict(
        list
    )  # So that we can at least know of collisions rather than silently dropping them.

    raw_input = open_raw_csv_ebird(csv_path)
    for line in raw_input:
        if line["category"] == "species":
            common_name = line["primary_com_name"]
            scientific_name = line["sci_name"]
            common_four_letter_code = name_to_4lc(common_name)
            scientific_four_letter_code = name_to_4lc(scientific_name)
            species_code = line["species_code"]
            short_codes = common_four_letter_code + scientific_four_letter_code
            taxon = Taxonomy(common_name, scientific_name, species_code, short_codes)
            common_map[common_name] = taxon
            scientific_map[scientific_name] = taxon
            code_map[species_code] = taxon
            if common_name == "Yellow-rumped Warbler":
                short_codes += ["MYWA", "AUWA"]
            for x in short_codes:
                short_map[x] += [taxon]
    return common_map, scientific_map, code_map, short_map


def common_name_to_banding(filename, include_non_sp=False):
    """
    Converts the banding code CSV file to a dictionary.
    Banding codes as per The Institute for Bird Populations downloaded from: http://www.birdpop.org/docs/misc/IBPAOU.zip
    Args:
        filename (str): csv file to open.
        include_non_sp (bool, optional): Whether ot not to include non-species taxa, such as subspecies or morphs. Defaults to False.
    Returns:
        dict: A dictionary of {"common name": "4 letter code"} pairs.
    """
    result = {}
    with open(filename, "r") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            if row["SP"] and not include_non_sp:
                continue
            common_name = row["COMMONNAME"]
            banding_code = row["SPEC"]
            result[common_name] = banding_code
    return result


def create_mapping_output(banding_csv_filename, ebird_csv_filename):
    """
    Simple script to open and parse the eBird teaxonomy csv and banding codes from Birdpop.
    Creates a json file called "short_codes.json" in the current directory.
    The result looks like {"Common Ostrich": {"common_name": "Common Ostrich", "scientific_name": "Struthio camelus", "species_code": "ostric2", "short_codes": ["COOS", "STCA"]}, ...}
    Args:
        banding_csv_filename (str): Path ot the banding code CSV.
        ebird_csv_filename (str): Path to the eBird taxonomy CSV.
    """
    # The banding codes don't always match up with eBird names. This dict maps the differences.
    # Format is "banding code name": "ebird name"
    # Right now eBird recognizes some of these and not others, or maps the same code to multiple species in the case of *certain* splits.
    # Anything that doesn't have a 1-1 mapping with ebird is currently ignored, but included for informational purposes.
    EBIRD_MAPPING = {
        "Fork-tailed Swift": "Pacific Swift",
        "Western Water-Rail": "Water Rail",
        # "Purple Swamphen": "???", Ignored as this has been split and is an invalid species.
        "Common Moorhen": "Eurasian Moorhen",
        "Western Marsh Harrier": "Eurasian Marsh-Harrier",
        "Checker-throated Antwren": "Checker-throated Stipplethroat",
        # "Paltry Tyrannulet": "???", Ignored. Another split that isn't a valid species in eBird.
        "Japanese Bush-Warbler": "Japanese Bush Warbler",
        "Japanese White-eye": "Warbling White-eye",  # eBird does not recognize this code.
        # "Hwamei": "???", Ignored. Another split that doesn't map to one species.
        "Stonechat": "European Stonechat",
        }

    common, _, _, _ = taxonomy_parse(ebird_csv_filename)

    banding_mapping = common_name_to_banding(banding_csv_filename)
    taxonmy_data = {k: v.as_dict() for k, v in common.items()}

    for name, banding_code in banding_mapping.items():
        name = EBIRD_MAPPING.get(name, name)
        try:
            taxonmy_data[name]["short_codes"] += [banding_code]
        except KeyError:
            pass  # Skip the ignored values.

    with open("short_codes.json", "w") as f:
        json.dump(taxonmy_data, f)


if __name__ == "__main__":
 create_mapping_output("list19p.csv", "eBird_Taxonomy_v2019.csv")
