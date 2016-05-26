import plinkio
import logging
import metax.Utilities
import metax.ThousandGenomesUtilities

def process_thousand_genomes(input_path, output_path, index, results, all_people, selected_people_by_id):
    names = metax.Utilities.hapNamesFromFolder(input_path)

    for name in names:
        logging.info("Processing %s", name)
        callback = PlinkBuilder(output_path, index, results, all_people, selected_people_by_id)
        loader = metax.ThousandGenomesUtilities.IMPUTELoader(input_path, name)
        loader.iterateOverFileDosage(callback)


class PlinkBuilder(object):
    def __init__(self, output_path, index, snps_by_gene, all_people, selected_people__by_id):
        self.output_path = output_path
        self.index = index
        self.snps_by_gene = snps_by_gene
        self.all_people = all_people
        self.selected_people_by_id = selected_people__by_id
        self.files = {}

        self.found_snps = 0
        self.snps_len = len(self.index)

    def __call__(self, hap_line, legend_line, row):
        rsid, valid, legend = metax.ThousandGenomesUtilities.checkLegend(legend_line, self.index)
        if not valid:
            return

        logging.log(3, "id %s found in white list", rsid)
        self.found_snps += 1
        percent = int(round(self.found_snps * 100.0 / self.snps_len))
        if percent > self.last_reported_percenteage:
            self.last_reported_percenteage = percent
            logging.log(9, "%d percent of snps found", percent)

        dosages = metax.ThousandGenomesUtilities.buildDosages(hap_line, self.all_people, self.selected_people_by_id)

        Now write to plink