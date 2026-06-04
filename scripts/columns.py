RECEPTOR = "Receptor"
CO_RECEPTOR = "Co-Receptor"
MOLECULE = "Molecule"
RESPONSE = "Response"
ASSAY = "Assay"
SOURCE = "Source"
BACKOFFICE = "Backoffice"

ORDER = "Order"
SPECIES = "Species"
GENE_NAME = "Gene Name"
UNIPROT_ID = "UniProt ID"
DATABASE = "Database"
IDENTITY = "Identity"  # New
MUTATION = "Mutation"
TAG = "Tag"
SEQUENCE = "Sequence"
SEQUENCE_REF = "Sequence_ref"

MOLECULE_NAME = "Molecule Name"
CID = "CID"
CAS = "CAS"
INCHIKEY = "InChIKey"
SMILES = "SMILES"
SOLVENT = "Solvent used for dilution"
MIXTURE = "Mixture"

RESPONSIVE = "Responsive"
PARAMETER = "Parameter"
VALUE = "Value"
UNIT = "Unit"
VALUE_NATURE = "Value Nature"
CONCENTRATION = "Concentration"
CONCENTRATION_UNIT = "Concentration Unit"
NBR_MEASUREMENTS = "Nbr. measurements"
MULTIMER = "Multimer responsive to control ligand (VUAA1 , …)"

EXPERIMENTAL_TECHNIQUE = "Experimental technique"
RECORDING_SYSTEM = "Recording System"
TYPE = "Type"
EXPRESSION_SYSTEM = "Expression system"
EXPRESSION_CELL_TYPE = "Expression cell type / line"
EXPRESSION_DRIVER = "Expression driver"
CO_TRANSFECTION = "Co-transfection"
ODOR_DELIVERY = "Odor delivery system"
STIMULATION_FLUX = "Stimulation Flux"
STIMULATION_FLUX_UNIT = "Stimulation Flux Unit"
MAIN_FLUX = "Main Flux"
MAIN_FLUX_UNIT = "Main Flux Unit"
STIMULATION_DURATION = "Stimulation duration"
STIMULATION_DURATION_UNIT = "Stimulation duration unit"
LARVA_OR_ADULT = "Larva or adult"

REFERENCE = "Reference"
DOI = "DOI"
REFERENCE_POSITION = "Reference Position"

COLLECTOR = "Collector"
EXTRA = "Extra"

_RECEPTOR_COLUMNS = [
    ORDER,
    SPECIES,
    GENE_NAME,
    UNIPROT_ID,
    MUTATION,
    TAG,
    SEQUENCE,
]

_CO_RECEPTOR_COLUMNS = [
    ORDER,
    SPECIES,
    GENE_NAME,
    UNIPROT_ID,
    MUTATION,
    SEQUENCE,
]

_RESPONSE_COLUMNS = [
    RESPONSIVE,
    PARAMETER,
    VALUE,
    UNIT,
    VALUE_NATURE,
    CONCENTRATION,
    CONCENTRATION_UNIT,
    NBR_MEASUREMENTS,
    MULTIMER,
]

_ASSAY_COLUMNS = [
    EXPERIMENTAL_TECHNIQUE,
    RECORDING_SYSTEM,
    TYPE,
    EXPRESSION_SYSTEM,
    EXPRESSION_CELL_TYPE,
    EXPRESSION_DRIVER,
    CO_TRANSFECTION,
    ODOR_DELIVERY,
    STIMULATION_FLUX,
    STIMULATION_FLUX_UNIT,
    MAIN_FLUX,
    MAIN_FLUX_UNIT,
    STIMULATION_DURATION,
    STIMULATION_DURATION_UNIT,
    LARVA_OR_ADULT,
]

COLUMNS_BY_GROUP: list[tuple[str, list[str]]] = [
    (RECEPTOR, _RECEPTOR_COLUMNS),
    (CO_RECEPTOR, _CO_RECEPTOR_COLUMNS),
    (MOLECULE, [MOLECULE_NAME, CID, CAS, INCHIKEY, SMILES, SOLVENT, MIXTURE]),
    (RESPONSE, _RESPONSE_COLUMNS),
    (ASSAY, _ASSAY_COLUMNS),
    (SOURCE, [REFERENCE, DOI, REFERENCE_POSITION]),
    (BACKOFFICE, [COLLECTOR, EXTRA]),
]

GROUPS_ORDER = [group for group, _ in COLUMNS_BY_GROUP]
