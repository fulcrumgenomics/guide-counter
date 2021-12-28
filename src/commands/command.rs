use anyhow::Result;
use enum_dispatch::enum_dispatch;

#[enum_dispatch]
pub trait Command {
    fn execute(&self) -> Result<()>;
}
