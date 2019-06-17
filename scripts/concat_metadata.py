import argparse
import pandas as pd

cast_to_int = ["residence_census_tract"]

region_mapping = {
    "connecticut": "northeast_usa",
    "maine": "northeast_usa",
    "massachusetts": "northeast_usa",
    "new_hampshire": "northeast_usa",
    "rhode_island": "northeast_usa",
    "vermont": "northeast_usa",
    "new_jersey": "northeast_usa",
    "new_york": "northeast_usa",
    "pennsylvania": "northeast_usa",
    "illinois": "midwest_usa",
    "indiana": "midwest_usa",
    "michigan": "midwest_usa",
    "ohio": "midwest_usa",
    "wisconsin": "midwest_usa",
    "iowa": "midwest_usa",
    "kansas": "midwest_usa",
    "minnesota": "midwest_usa",
    "missouri": "midwest_usa",
    "nebraska": "midwest_usa",
    "north_dakota": "midwest_usa",
    "south_dakota": "midwest_usa",
    "delaware": "south_usa",
    "florida": "south_usa",
    "georgia": "south_usa",
    "maryland": "south_usa",
    "north_carolina": "south_usa",
    "south_carolina": "south_usa",
    "virginia": "south_usa",
    "district_of_columbia": "south_usa",
    "west_virginia": "south_usa",
    "alabama": "south_usa",
    "kentucky": "south_usa",
    "mississippi": "south_usa",
    "tennessee": "south_usa",
    "arkansas": "south_usa",
    "louisiana": "south_usa",
    "oklahoma": "south_usa",
    "texas": "south_usa",
    "arizona": "west_usa",
    "colorado": "west_usa",
    "idaho": "west_usa",
    "montana": "west_usa",
    "nevada": "west_usa",
    "new_mexico": "west_usa",
    "utah": "west_usa",
    "wyoming": "west_usa",
    "alaska": "west_usa",
    "california": "west_usa",
    "hawaii": "west_usa",
    "oregon": "west_usa",
    "washington": "west_usa",
    "canada": "canada",
    "anguilla": "central_america",
    "bahamas": "central_america",
    "barbados": "central_america",
    "belize": "central_america",
    "bermuda": "central_america",
    "cayman_islands": "central_america",
    "costa_rica": "central_america",
    "cuba": "central_america",
    "dominica": "central_america",
    "dominican_republic": "central_america",
    "el_salvador": "central_america",
    "grenada": "central_america",
    "guadeloupe": "central_america",
    "guatemala": "central_america",
    "haiti": "central_america",
    "honduras": "central_america",
    "jamaica": "central_america",
    "martinique": "central_america",
    "mexico": "central_america",
    "nicaragua": "central_america",
    "panama": "central_america",
    "puerto_rico": "central_america",
    "saint_lucia": "central_america",
    "saint_vincent_and_the_grenadines": "central_america",
    "trinidad_and_tobago": "central_america",
    "turks_and_caicos": "central_america",
    "usa": "?"
}

def select(file, mergeby, fields):
    entries = pd.read_csv(file, sep='\t')
    mapping = {}
    for index, row in entries.iterrows():
        if row[mergeby]:
            values = []
            for field in fields:
                if field in row:
                    if str(row[field]) != 'nan':
                        if field == "region":
                            if "country" in row and "division" in row:
                                if row["division"] in region_mapping:
                                    values.append(region_mapping[row["division"]])
                                elif row["country"] in region_mapping:
                                    values.append(region_mapping[row["country"]])
                                else:
                                    values.append(row["region"])
                            else:
                                values.append(row["region"])
                        else:
                            if field in cast_to_int:
                                values.append(str(int(row[field])))
                            else:
                                values.append(str(row[field]))
                    else:
                        values.append("?")
                else:
                    values.append("?")
            mapping[row[mergeby]] = values
    for key, values in mapping.items():
        print(str(key) + "\t" + "\t".join(values))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Concatenate multiple tsvs, merging specified columns",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--files', nargs='+', type=str, required=True, help="list of tsv files")
    parser.add_argument('--mergeby', type=str, help="column name to merge on")
    parser.add_argument('--fields', nargs='+', type=str, help="column names to include")
    args = parser.parse_args()

    print(args.mergeby + "\t" + "\t".join(args.fields))
    for file in args.files:
        select(file, args.mergeby, args.fields)
