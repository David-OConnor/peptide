use core::f64::consts::TAU;

use egui;

const WINDOW_TITLE: &str = "Peptide info";
const WINDOW_SIZE_X: f32 = 900.0;
const WINDOW_SIZE_Y: f32 = 600.0;

// Note: This is draggable.
const SIDE_PANEL_SIZE: f32 = 400.;

// todo: Quick and dirty here, just like in `render_gpu`. Probably avoidable
// todo by using Fn traits instead of `fn` pointers.
pub static mut PROT_NAME: &'static str = "";
pub static mut PDB_IDENT: &'static str = "";
pub static mut ACTIVE_RES_ID: usize = 1;
pub static mut ACTIVE_RES_AA_NAME: &'static str = "";
pub static mut ACTIVE_RES_PHI: f64 = 0.;
pub static mut ACTIVE_RES_PSI: f64 = 0.;
pub static mut ACTIVE_RES_OMEGA: f64 = 0.;

/// This function draws the (immediate-mode) GUI. We're currently editing it
/// from the main program by modifying the `static mut` variables above.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html#method.heading)
pub fn draw_ui(ctx: &egui::Context) {
    unsafe {
        let panel = egui::SidePanel::left(0)  // ID must be unique among panels.
            .default_width(SIDE_PANEL_SIZE);

        panel.show(ctx, |ui| {
            // println!("{:?}", ui.spacing());
            ui.spacing_mut().item_spacing = egui::vec2(10.0, 12.0);

            // ui.label("Protein: ".to_owned().push_str(prot_name));
            ui.heading(format!("Protein: {PROT_NAME}. PDB: {PDB_IDENT}"));

            ui.label(format!("Active Residue: {ACTIVE_RES_ID}"));

            ui.label(ACTIVE_RES_AA_NAME);

            ui.horizontal(|ui| {
                // ui.text_edit_singleline(&mut aa_name);
            });

            // todo: Put the display in terms of Tau.
            ui.add(egui::Slider::new(&mut ACTIVE_RES_PHI, 0.0..=TAU).text("φ"));
            ui.add(egui::Slider::new(&mut ACTIVE_RES_PSI, 0.0..=TAU).text("ψ"));
            ui.add(egui::Slider::new(&mut ACTIVE_RES_OMEGA, 0.0..=TAU).text("ω"));

            //
            // if ui.button("Click each year").clicked() {
            //     // Perform action here.
            // }
        });
    }
}
