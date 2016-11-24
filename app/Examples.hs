import           Graphics.Vty

main :: IO ()
main = do
    cfg <- standardIOConfig
    vty <- mkVty cfg
    let line0 = string (defAttr ` withForeColor ` green) "first line"
        line1 = string (defAttr ` withBackColor ` blue) "second line"
        img = line0 <-> line1
        pic = picForImage img
    update vty pic
    e <- nextEvent vty
    shutdown vty
    print ("Last event was: " ++ show e)

-- main = putStrLn "My name is William Rowan Hamilton!"
