use winit::{
    event::{DeviceEvent, Event, KeyboardInput, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    window::{Window, WindowBuilder},
};

const WINDOW_TITLE: &str = "Peptide info";
const WINDOW_SIZE_X: f32 = 900.0;
const WINDOW_SIZE_Y: f32 = 600.0;

pub fn setup() {
    let event_loop = EventLoop::new();
    let window = WindowBuilder::new()
        .with_title(WINDOW_TITLE)
        .with_inner_size(winit::dpi::LogicalSize::new(WINDOW_SIZE_X, WINDOW_SIZE_Y))
        .build(&event_loop)
        .unwrap();

    event_loop.run(move |event, _, control_flow| {
        *control_flow = ControlFlow::Poll;

        match event {
            Event::MainEventsCleared => window.request_redraw(),
            Event::DeviceEvent { event, .. } => {}
            Event::WindowEvent {
                ref event,
                window_id,
                // } if window_id == window.id() && !state.input(event) => {
            } if window_id == window.id() => {
                match event {
                    // todo: Put back for window-closing.
                    // #[cfg(not(target_arch="wasm32"))]
                    WindowEvent::CloseRequested
                    // | WindowEvent::KeyboardInput {
                    //     input:
                    //     KeyboardInput {
                    //         state: ElementState::Pressed,
                    //         virtual_keycode: Some(VirtualKeyCode::Escape),
                    //         ..
                    //     },
                    //     ..
                    => *control_flow = ControlFlow::Exit,
                    WindowEvent::Resized(physical_size) => {
                        // state.resize(*physical_size);
                    }
                    WindowEvent::ScaleFactorChanged { new_inner_size, .. } => {
                        // state.resize(**new_inner_size);
                    }
                    _ => {}
                }
            }

            Event::RedrawRequested(window_id) if window_id == window.id() => {}
            _ => {}
        }
    });
}
