use egui;

const WINDOW_TITLE: &str = "Peptide info";
const WINDOW_SIZE_X: f32 = 900.0;
const WINDOW_SIZE_Y: f32 = 600.0;

// todo: Quick and dirty here, just like in `render_gpu`. Probably avoidable
// todo by using Fn traits instead of `fn` pointers.
pub static mut PROT_NAME: &'static str = "";
pub static mut PDB_IDENT: &'static str = "";
pub static mut ACTIVE_RES_ID: usize = 1;
pub static mut ACTIVE_RES_AA_NAME: &'static str = "";

/// This function draws the (immediate-mode) GUI. We're currently editing it
/// from the main program by modifying the `static mut` variables above.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html#method.heading)
pub fn draw_ui(ctx: &egui::Context) {
    unsafe {
        // let mut aa_name = "Amino Acid: ".to_owned().push_str("");

        let mut age = 20;

        let panel = egui::SidePanel::left(0); // ID must be unique among panels.
        panel.show(ctx, |ui| {
            // ui.heading("Peptide");

            // ui.label("Protein: ".to_owned().push_str(prot_name));
            ui.label(format!("Protein: {PROT_NAME}. PDB: {PDB_IDENT}"));

            ui.label(format!("Active Residue: {ACTIVE_RES_ID}"));

            ui.label(ACTIVE_RES_AA_NAME);

            ui.horizontal(|ui| {
                ui.label("Label with edit?");
                // ui.text_edit_singleline(&mut aa_name);
            });

            ui.add(egui::Slider::new(&mut age, 0..=120).text("age"));

            if ui.button("Click each year").clicked() {
                // Perform action here.
            }
        });
    }
}
