use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about = "kv-mer", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// Sketch the given sequencing files into kv-mer sketches.
    //#[clap(display_order = 1)]
    //Sketch(AnalyzeArgs),

    /// Analyze a given sequencing file.
    #[clap(display_order = 1)]
    Analyze(AnalyzeArgs),
}

#[derive(Args, Default)]
pub struct SketchArgs {
    #[clap(multiple=true, help_heading = "INPUT", help = "fasta/fastq files; gzip optional.")]
    pub files: Vec<String>,

    #[clap(short, default_value_t = 21, help_heading = "ALGORITHM", help ="Length of keys.")]
    pub k: u8,

    #[clap(short, default_value_t = 6, help_heading = "ALGORITHM", help ="Length of values.")]
    pub v: u8,

    #[clap(short, default_value_t = 200, help_heading = "ALGORITHM", help = "Subsampling rate.")]
    pub c: usize,

    
}


#[derive(Args, Default)]
pub struct AnalyzeArgs {
    #[clap(multiple=true, help_heading = "INPUT", help = "fasta/fastq files; gzip optional.")]
    pub files: Vec<String>,

    #[clap(short, default_value_t = 21, help_heading = "ALGORITHM", help ="Length of keys.")]
    pub k: u8,

    #[clap(short, default_value_t = 6, help_heading = "ALGORITHM", help ="Length of values.")]
    pub v: u8,

    #[clap(short, default_value_t = 200, help_heading = "ALGORITHM", help = "Subsampling rate.")]
    pub c: usize,

    #[clap(short, default_value_t = 200, help_heading = "ALGORITHM", help = "Threshold for consensus.")]
    pub threshold: u32,
}