mod foreign_call;
mod handlers;
mod ops;

use jsonrpsee::server::{RpcModule, Server};
use std::net::SocketAddr;
use tracing_subscriber::util::SubscriberInitExt;

use serde::Deserialize;
use serde_json::{json, Value};

use crate::foreign_call::ForeignCallParam;
use crate::handlers::{
    handle_get_pairing_witnesses, handle_is_third_root, handle_random_third_root,
    handle_third_root, handle_witness_gen,
};

// SPIN UP THE SERVER
#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let filter = tracing_subscriber::EnvFilter::try_from_default_env()?
        .add_directive("jsonrpsee[method_call{name = \"say_hello\"}]=trace".parse()?);
    tracing_subscriber::FmtSubscriber::builder()
        .with_env_filter(filter)
        .finish()
        .try_init()?;

    let _server_addr = run_server().await?;

    Ok(())
}

fn print_type<T>(_: &T) {
    println!("{:?}", std::any::type_name::<T>());
}

#[derive(Debug, Deserialize)]
struct RequestData {
    session_id: u64,
    function: String,
    inputs: Vec<ForeignCallParam<String>>,
    root_path: String,
    package_name: String,
}

#[derive(Debug, Deserialize)]
struct Requests(Vec<RequestData>); // Wrap it in a struct to handle the array

pub(crate) fn handle_unknown_function(_input: &RequestData) -> Value {
    println!("oops");
    json!(vec![String::from("oops")])
}

async fn run_server() -> anyhow::Result<SocketAddr> {
    let server = Server::builder()
        .build("127.0.0.1:3000".parse::<SocketAddr>()?)
        .await?;
    let mut module = RpcModule::new(());

    module.register_method("say_hello", |_, _, _| "hello, world")?;

    module.register_method("resolve_foreign_call", |params, _, _| {
        // println!("\n\nNEW REQUEST!!!");
        // print_type(&params);
        // println!("params{:?}", params);

        let response: Value = if let Some(json_string) = params.as_str() {
            // Deserialize the JSON string into the Requests struct:
            let requests: Requests =
                serde_json::from_str(&json_string).expect("Failed to parse JSON");

            let request = &requests.0[0];

            // println!("Request function{:?}", request.function);

            let result: Value = match request.function.as_str() {
                "witness_gen" => handle_witness_gen(&request.inputs),
                "third_root" => handle_third_root(&request.inputs),
                "is_third_root" => handle_is_third_root(&request.inputs),
                "random_third_root" => handle_random_third_root(&request.inputs),
                "get_pairing_witnesses" => handle_get_pairing_witnesses(&request.inputs),
                _ => handle_unknown_function(&request),
            };

            result
        } else {
            println!("No parameters provided");
            json!(vec![String::from("Bad query")])
        };

        response
    })?;

    let addr = server.local_addr()?;
    let handle = server.start(module);

    println!("Server is running on 127.0.0.1:3000");

    // In this example we don't care about doing shutdown so let's it run forever.
    // You may use the `ServerHandle` to shut it down or manage it yourself.
    // tokio::spawn(handle.stopped());

    // Keep the server running until it's interrupted
    handle.stopped().await;

    Ok(addr)
}
