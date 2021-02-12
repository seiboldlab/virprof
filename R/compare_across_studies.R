#!/usr/bin/env Rscript
library(here)
library(tidyverse)
library(janitor)
library(patchwork)
library(ggrepel)
library(ggforce)
library(ggbeeswarm)
library(cowplot)
library(ggpubr)
library(GGally)
library(magrittr)


library(httr)
library(xml2)
library(progress)

NCBI_genome_size_check <- function(taxids) {
    API_SERVER <- "https://api.ncbi.nlm.nih.gov"
    API_PATH <- "genome/v0/expected_genome_size"
    unique_taxids <- unique(taxids)
    pb <- progress::progress_bar$new(total=length(unique_taxids))
    pb$tick(0)
    check_one <- function(taxid) {
        url <- httr::modify_url(API_SERVER, path=API_PATH,
                                query = list(species_taxid=taxid))
        response <- GET(url)
        if (httr::status_code(response) != 200) {
            stop("status code")
        }
        pb$tick()
        text = httr::content(response, "text", encoding = "utf-8")
        dx <- xml2::read_xml(text)
        xml2::xml_ns_strip(dx)
        nodes <- xml2::xml_find_all(dx, "/genome_size_response/*[not(self::input)]")
        res <- set_names(
            sapply(nodes, xml2::xml_text),
            sapply(nodes, xml2::xml_name)
        )
    }
    unique_results <- lapply(unique_taxids, check_one)
    results <- unique_results[match(taxids, unique_taxids)]
    do.call(bind_rows, results) %>%
        as_tibble() %>%
        type_convert()
}

load_data <- function(.data, .f, .file_pattern="%s", .id=NULL, ...) {
    set_names(
        sprintf(.file_pattern, .data),
        .data
    ) %>%
        map_df(.f, .id=.id, ...)
}

studies <- c("gala", "daps", "inspire", "enigma")

genes_of_interest <- c(
    "CCL8", "CXCL11"
)


## Load Meta Data
units <- load_data(studies, read_csv, here("projects/%s.csv"), .id="study") %>%
    mutate(sample=if_else(is.na(sample), unit, sample)) %>%
    mutate_at(vars(study, unit, sample), as.factor)

##units %>% print(n=50900)

## Load Gene Expression Data
expression_data <- load_data(
    studies,
    function(x) {
        readRDS(x)$dds_norm %>%
                     rownames_to_column("gene") %>%
                     pivot_longer(-gene, names_to="sample")
    },
    here("expression_data", "%s.rds"),
    .id="study")


## Load Virus Detection Data
rviruses_cols <- cols(
    unit = col_character(),
    sample = col_character(),
    words = col_character(),
    taxname = col_character(),
    sacc = col_character(),
    stitle = col_character(),
    staxids = col_character(),
    log_evalue = col_double(),
    pident = col_double(),
    qlen = col_double(),
    alen = col_double(),
    n_frag = col_double(),
    qaccs = col_character(),
    qranges = col_character(),
    sranges = col_character(),
    numreads = col_double(),
    taxid = col_double(),
    lineage = col_character(),
    host = col_character()
)

data_path <- paste(collapse="", c("%s",
                                  ".remove_human_reads",
                                  ".assemble_and_map",
                                  ".ref_hg38.annotate_blastE2MegaBest.",
                                  "extract_seqsNomatch.ref_NT.annotate_blastE10Best.",
                                  "ref_NcbiTaxonomy.classify"
                                  ))
filenames <- here(data_path, "result_resp_viruses.csv"),

rviruses <- load_data(
    studies, read_csv,
    here("results/%s_resp_viruses.csv"),
    col_types=rviruses_cols,
    .id="study") %>%
    mutate_at(vars(study, unit, sample, lineage, host),
              as.factor)

## Determine Virus Genome Sizes
genome_sizes <- NCBI_genome_size_check(rviruses$taxid)

manual_sizes <- tribble(
    ~ organism_name,  ~ species_taxid, ~ len, ~acc,
    "Human respirovirus 1",     12730, 15600, "NC_003461",
    "Human coronavirus NL63",  277944, 27553, "NC_005831",
    "Betacoronavirus 1",       694003, 30741, "NC_006213", ## OC43
    "Human respirovirus 3",     11216, 15462, "NC_001796", # parainfluenza
    "Human coronavirus 229E",   11137, 27317, "NC_002645",
    "Human coronavirus HKU1",  290028, 29926, "NC_006577",
    "Enterovirus D",           138951,  7390, "NC_001430", #D70
    "Human rhinovirus AMS323", 433730,  7160, "EF155421",
    "Influenza B virus",        11520, 2313 + 2204 + 1882 + 1841 + 1557 +
                                       1191 + 1096 + 2368, "NC_002204", 
    "Human rhinovirus sp.",    169066,  7208, "X01087.1" #HRV-14
)

genome_sizes %<>%
    left_join(select(manual_sizes, species_taxid, len)) %>%
    mutate(
        expected_ungapped_length = if_else(
            is.na(expected_ungapped_length),
            len,
            expected_ungapped_length
        )
    ) %>%
    select(-len)

if (any(is.na(genome_sizes$expected_ungapped_length))) {
    genome_sizes[is.na(genome_sizes$expected_ungapped_length), ]
    stop("Missing sizes for genomes")
}

rviruses <- rviruses %>%
    bind_cols(
        select(genome_sizes,
               species_taxid, expected_ungapped_length)
    ) %>%
    rename(genome_length = expected_ungapped_length)

## Load raw read counts
raw_read_counts <- paste0(
    "grep 'Input:' ",
    paste0(collapse=" ", here(""), studies, ".trim_bbmapAQ10/*.log"),
    "| sed 's/.*\\/\\(.*\\)\\.trim.*\\/\\(.*\\).log:\\w*\\W*\\([0-9]*\\) reads.*/\\1,\\2,\\3/'"
) %>%
    system(intern = TRUE) %>%
    read_csv(col_names = c("study", "unit", "raw_reads"))

## Load filtered read counts
filtered_read_counts <- paste0(
    "grep 'Result:' ",
    paste0(collapse=" ", here(""), studies, "*.dust_bbmapE??/*.log"),
    "| sed 's/.*\\/\\([^_]*\\)\\..*dust.*\\/\\(.*\\).log:\\w*\\W*\\([0-9]*\\) reads.*/\\1,\\2,\\3/'"
) %>%
    system(intern = TRUE) %>%
    read_csv(col_names = c("study", "unit", "filtered_reads"))

## Fixup sample names
expression_data %<>%
    mutate(
        sample = sub("DAPS-", "DAPS", sample)
    )

renamed_samples = c(
    HR5579 = "HR5519",
    SJ1241 = "HR1241",
    HR1212 = "HR5125",
    HR1250 = "HR1251",
    SJ5398 = "HR5398",
    HR1206 = "HR1206a"
)

units %<>%
    mutate(
        sample = sub("DAPS([^_]*)_.*", "DAPS\\1", sample),
        sample = as.factor(sample),
        sample = fct_recode(sample, !!!renamed_samples)
    )

rviruses %<>%
    mutate(
        sample = sub("DAPS([^_]*)_.*", "DAPS\\1", sample),
        sample = as.factor(sample),
        sample = fct_recode(sample, !!!renamed_samples),
        unit = sub(".virus.csv$", "", unit) # enigma
    )



## Match sample names

matched <- expression_data %>%
    group_by(study, sample) %>%
    summarize(have_expr=TRUE) %>%
    ungroup() %>%
    mutate_if(is.character, as.factor) %>%
    full_join(select(units, study, sample, unit)) %>%
    mutate(
        have_meta = !is.na(unit),
        have_expr = !is.na(have_expr)
    ) %>%
    full_join(select(rviruses, study, sample, unit, n_frag)) %>%
    mutate(
        have_virus = !is.na(n_frag)
    ) %>%
    select(-unit, -n_frag) %>%
    filter(!have_meta | !have_expr)
matched %>% print(n=500)

data <- units %>%
#    full_join(raw_read_counts, by=c("study", "unit")) %>%
#    full_join(filtered_read_counts, by=c("study", "unit")) %>%
    left_join(rviruses, by=c("study", "sample", "unit")) %>%
    mutate(
        Virus = taxname,
        Virus = sub(" virus.*", "", Virus),
        Virus = sub("rhinovirus sp.", "rhinovirus B", Virus), # CAUTION: this is specific to one sample
        Virus = sub("rhinovirus 1B", "rhinovirus A", Virus), # 1B is an A, really
        Virus = sub("(Human )?[Rr]hinovirus ([A-Z]).*", "HRV \\2", Virus),
        Virus = sub("Human orthopneumovirus", "RSV", Virus),
        Virus = sub("Human metapneumovirus", "HMPV", Virus),
        Virus = sub("Human respirovirus", "HPIV", Virus),
        Virus = sub("Human coronavirus", "CoV", Virus),
        Virus = sub("(Human )?[Ee]nterovirus", "EV", Virus),
        study = if_else(study == "enigma",
                        gsub(" .*_", " ", paste(study, sample)),
                        as.character(study)),
        ) %>%
    separate(Virus, c("Virus", "Strain"), fill="right")



#' Graphs showing:
#'
#' Breakdown Positive vs Negative
g1 <- data %>%
    group_by(study, sample) %>%
    summarize(
        Infection = if(all(is.na(taxname))) "negative" else "positive"
    ) %>%
    group_by(study, Infection) %>%
    summarize(n=n())  %>%
    group_by(study) %>%
    mutate(p=n/sum(n)) %>%
    ggplot(aes(x=1, y=p, fill=Infection, label=p)) +
    facet_grid(~ study) +
    geom_bar(
        width=1,
        color="white",
        stat="identity"
    ) +
    geom_text(
        aes(label = sprintf("%02.1f%%\n(%i)", p*100, n)),
        stat="identity",
        position=position_fill(vjust=.5)
    ) +
    theme_minimal() +
    theme(
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid   = element_blank(),
        axis.ticks   = element_blank()
    ) +
    coord_polar(
        theta = "y",
        start = 0
    ) + 
    ggtitle("Viral Infection Rate")

g2 <- data %>%
    filter(!is.na(Virus)) %>%
    group_by(study, sample) %>%
    summarize(Virus = paste(collapse="+", unique(Virus))) %>%
    group_by(study, Virus) %>%
    summarize(n = n()) %>%
    mutate(p = n/sum(n)) %>%
    ggplot(aes(x=Virus, fill=Virus, y=p)) +
    geom_col(
       position = position_dodge2(preserve = "single")
    ) +
    geom_text_repel(
        aes(label = sprintf("%02.1f%%\n(%i)", p*100, n)),
        size=10/.pt,
        direction = "y",
        box.padding = 0.1
    ) +
    facet_grid(~ study, space = "free", scales = "free") +
    theme(
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5)
    ) +
    guides(
        fill=FALSE
    ) +
    scale_y_continuous(
        labels = scales::label_percent()
    ) + ggtitle("Virus Distribution")

g3 <- data %>%
    filter(Virus %in% c("HRV")) %>%
    group_by(study, sample) %>%
    summarize(Strain = paste(collapse="+", unique(Strain))) %>%
    group_by(study, Strain) %>%
    summarize(n = n()) %>%
    group_by(study) %>%
    mutate(p = n/sum(n)) %>%
    ggplot(aes(x=Strain, fill=Strain, y=p)) +
    geom_col(
        position = position_dodge2(preserve = "single")
    ) + 
    geom_text_repel(
        aes(label = sprintf("%02.1f%%\n(%i)", p*100, n)),
        direction="y",
        size=10/.pt
    ) +
    facet_grid(~ study, space = "free_x", scales = "free", drop=TRUE) +
    theme(
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5)
    ) +
    guides(
        fill=FALSE
    ) +
    scale_y_continuous(
        labels = scales::label_percent()
    )  + ggtitle("Rhinovirus Distribution")

g4 <- data %>%
    filter(Virus %in% c("CoV")) %>%
    group_by(study, sample) %>%
    summarize(Strain = paste(collapse="+", unique(Strain))) %>%
    group_by(study, Strain) %>%
    summarize(n = n()) %>%
    group_by(study) %>%
    mutate(p = n/sum(n)) %>%
    ggplot(aes(x=Strain, fill=Strain, y=p)) +
    geom_col(
        position = position_dodge2(preserve = "single")
    ) + 
    geom_text_repel(
        aes(label = sprintf("%02.1f%%\n(%i)", p*100, n)),
        direction="y",
        size=10/.pt
    ) +
    facet_grid(~ study, space = "free_x", scales = "free", drop=TRUE) +
    theme(
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5)
    ) +
    guides(
        fill=FALSE
    ) +
    scale_y_continuous(
        labels = scales::label_percent()
    )  + ggtitle("Coronavirus Distribution")
ggsave(
    here("FIGURES", "infection_breakdowns.pdf"),
    (g1/g2/g3/g4 & theme(axis.title.x = element_blank(),
                         axis.title.y = element_blank())) + 
    plot_layout(widths=1),
    width=11, height=11
)



#' 2.Generate
#' - the viral reads,
#' - viral read counts normalized for total RNA-seq reads and viral genome size
#' - viral read counts normalized for total unmapped reads (going in to the contig assembly step) and
#' - viral genome size for all infected subjects

## - Gather per-sample counts:
##   - raw reads
##   - non-human reads
## - Determine Genome sizes
##

g5 <- data %>%
    filter(!is.na(Virus)) %>%
    group_by(study, sample) %>%
    summarize(
        Virus = paste0(collapse="+", unique(Virus)),
        Total = sum(raw_reads),
        Nonhuman = sum(filtered_reads),
        Viral = sum(numreads),
        PctViral = Viral / Total * 100,
        PctViralNonhuman = Viral / Nonhuman * 100
    ) %>%
    pivot_longer(
        c(Total, Nonhuman, Viral, PctViral, PctViralNonhuman),
        names_to = "norm"
    ) %>%
    mutate(
        norm = fct_relevel(norm,
                           "Total", "Nonhuman", "Viral",
                           "PctViral", "PctViralNonhuman")
    ) %>%
    ggplot(aes(
        color=Virus, 
        x=study, y=value
    )) +
    facet_wrap(~ norm, scales = "free", drop=TRUE) +
    geom_boxplot(color="grey", outlier.size = 0) +
    geom_quasirandom() +
    scale_y_log10() + 
    theme_cowplot()


log10_points <- function(data, mapping, ...) {
    ggally_autopoint(data, mapping, ...) + scale_x_log10() + scale_y_log10()
}
log10_diagonal <- function(data, mapping, ...) {
    ggally_densityDiag(data, mapping, ...) + scale_x_log10()
}

    
## get genome sizes and normalize

g6 <- data %>%
    filter(!is.na(Virus)) %>%
    group_by(study, sample) %>%
    summarize(
        Virus = paste0(collapse="+", unique(Virus)),
        Total = sum(raw_reads),
        Nonhuman = sum(filtered_reads),
        Genome = max(genome_length),
        Viral = sum(numreads),
        PctViral = Viral / Total * 100,
        PctViralNonhuman = Viral / Nonhuman * 100,
        ViralPerGenome = Viral / Genome * 100,
        RPKM_All = Viral / (Genome/1000) / (Total/1000000),
        RPKM_Nonhuman = Viral / (Genome/1000) / (Nonhuman/1000000),
        ) %>%
    ggpairs(
        columns = c("Viral", "PctViral", #"PctViralNonhuman",
                    "RPKM_All", "RPKM_Nonhuman"),
        mapping = aes(color=Virus),
        lower = list(
            continuous = log10_points
        ),
        diag = list(
            continuous = wrap(log10_diagonal, alpha=.5)
        ),
        upper = list(
            continuous = "blank"
        ),
        legend=5
    )
ggsave(here("FIGURES", "scaling.pdf"), g6, width=12, height=9)


log10_scatter <- function(data, mapping, ...) {
    ggscatter(
        data = data,
        x = "mapping$x"
        y = "mapping$y"
    )
}
 
g7 <- data %>%
    filter(!is.na(Virus)) %>%
    group_by(study, sample) %>%
    summarize(
        Virus = paste0(collapse="+", unique(Virus)),
        Total = sum(raw_reads),
        Nonhuman = sum(filtered_reads),
        Genome = max(genome_length),
        Viral = sum(numreads),
        PctViral = Viral / Total * 100,
        PctViralNonhuman = Viral / Nonhuman * 100,
        ViralPerGenome = Viral / Genome * 100,
        RPKM_All = Viral / (Genome/1000) / (Total/1000000),
        RPKM_Nonhuman = Viral / (Genome/1000) / (Nonhuman/1000000),
        ) %>%
    left_join(
        filter(expression_data, gene %in% genes_of_interest)
    ) %>%
    filter(sample != "C5046X4") %>%
    pivot_wider(names_from = "gene") %>%
    ggduo(
        columnsX = genes_of_interest,
        columnsY = c("Viral", "PctViral", "RPKM_All", "RPKM_Nonhuman"),
        mapping = aes(color=Virus, shape=study),
        legend=5,
        types = list(
            continuous = log10_scatter
        )
    )
ggsave(here("FIGURES", "tmp.pdf"), g7, width=12, height=9)


g7 <- data %>%
    filter(!is.na(Virus)) %>%
    group_by(study, sample) %>%
    summarize(
        Virus = paste0(collapse="+", unique(Virus)),
        raw_reads = sum(raw_reads),
        filtered_reads = sum(filtered_reads),
        genome_length = max(genome_length),
        `Read Count` = sum(numreads),
        `% Reads` = sum(numreads) / raw_reads * 100,
        `% Filt. Reads` = sum(numreads) / filtered_reads * 100,
        `Reads / Genome` = sum(numreads) / genome_length,
        `% Reads / Genome` = sum(numreads) / raw_reads * 100 / genome_length,
        `% Filt. Reads / Genome` = sum(numreads) / filtered_reads * 100 / genome_length,
        ) %>%
    left_join(
        filter(expression_data, gene %in% genes_of_interest)
    ) %>%
    filter(sample != "C5046X4") %>%
    pivot_longer(cols = c(
                      `Read Count`,
                     `% Reads`,
                     `% Filt. Reads`,
                     `Reads / Genome`,
                     `% Reads / Genome`,
                     `% Filt. Reads / Genome`
                 ),
                 values_to = "Viral Load",
                 names_to = "count_type"
                 ) %>%
    rename(
        Expression = value
    ) %>%
    ggscatter(
        y = "Expression",
        color = "Virus",
        shape = "study",
        xscale = "log10",
        x= "Viral Load",
        yscale = "log10",
        facet.by = c("gene + Virus", "count_type"),
        add = "reg.line",
        #add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
        conf.int = TRUE, 
        cor.coef = TRUE, 
#        cor.coeff.args = list(
#            method = "pearson",
#            label.x = 3,
#            label.sep = "\n"
#        ),
        scales="free"
    )
ggsave(here("FIGURES", "tmp.pdf"), g7, width=16, height=27)






ggpairs(data, columns = (")


#' 3. Scatter plot all infected subjects’ viral reads numbers (using
#' the 3 methods) versus gene expression of CXCL11 and CCL8. Make
#' separate plots for each cohort.
